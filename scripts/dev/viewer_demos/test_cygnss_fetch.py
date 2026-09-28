"""Offline data-access contract checks; no credentials or flight data required."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

spec = importlib.util.spec_from_file_location("fetch", Path(__file__).with_name("fetch_cygnss_l1_states.py"))
f = importlib.util.module_from_spec(spec); spec.loader.exec_module(f)

def inventory():
    return [dict(title=f"cyg{fm:02d}.ddmi.s202506{d:02d}-000000-test", opendap="https://opendap.earthdata.nasa.gov/test") for fm in f.EXPECTED_SPACECRAFT for d in range(6,10)]

class AccessTests(unittest.TestCase):
    def test_complete_inventory(self):
        self.assertEqual(f.validate_inventory(inventory()),[1,2,3,4,5,7,8])
    def test_missing_craft_is_not_success(self):
        with self.assertRaises(RuntimeError): f.validate_inventory(inventory()[:-4])
    def test_duplicate_granule_rejected(self):
        with self.assertRaises(RuntimeError): f.validate_inventory(inventory()+inventory()[:1])
    def test_untrusted_metadata_url_rejected(self):
        g=inventory();g[0]['opendap']='https://example.org/test'
        with self.assertRaises(RuntimeError): f.validate_inventory(g)
    def test_credential_domain(self):
        for url in ['http://opendap.earthdata.nasa.gov/x','https://earthdata.nasa.gov.evil.test/x','https://evil.test/x']:
            with self.assertRaises(RuntimeError): f._get(url,{'Authorization':'Bearer synthetic'})
        self.assertTrue(f.earthdata_url('https://opendap.earthdata.nasa.gov/x'))
    def test_redirect_guard(self):
        req=f.urllib.request.Request('https://opendap.earthdata.nasa.gov/x',headers={'Authorization':'Bearer synthetic'})
        with self.assertRaises(RuntimeError): f.EarthdataRedirect().redirect_request(req,None,302,'',{},'https://example.org/x')
    def test_missing_token_message(self):
        with tempfile.TemporaryDirectory() as d, patch.object(f,'TOKEN_PATH',Path(d)/'absent'):
            with self.assertRaisesRegex(RuntimeError,'Earthdata token unavailable'): f._bearer_header()
    def test_truncated_dap(self):
        for blob in [b'',b'\x01',b'\x05\x00\x00\x05ab']:
            with self.assertRaises(RuntimeError): f.dap4_dechunk(blob)
    def test_bad_cache_is_retried(self):
        with tempfile.TemporaryDirectory() as d:
            cache=Path(d);g={'title':'cyg01.test','opendap':'https://opendap.earthdata.nasa.gov/x'}
            (cache/'cyg01.test.dap').write_bytes(b'partial')
            with patch.object(f,'granule_rows',side_effect=[RuntimeError('bad'),{}]), patch.object(f,'_get',return_value=b'good') as get, patch.object(f,'_bearer_header',return_value={}):
                blob,n=f.fetch_granule_variables(g,cache,0)
            self.assertEqual((blob,n),(b'good',4));self.assertEqual(get.call_count,1)
            self.assertEqual((cache/'cyg01.test.dap').read_bytes(),b'good')
    def test_bad_response_not_cached(self):
        with tempfile.TemporaryDirectory() as d:
            with patch.object(f,'granule_rows',side_effect=RuntimeError('bad')), patch.object(f,'_get',return_value=b'bad'), patch.object(f,'_bearer_header',return_value={}):
                with self.assertRaises(RuntimeError): f.fetch_granule_variables({'title':'cyg01.test','opendap':'https://opendap.earthdata.nasa.gov/x'},Path(d),0)
            self.assertEqual(list(Path(d).iterdir()),[])
    def test_atomic_replacement(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a';p.write_bytes(b'old');f.atomic_bytes(p,b'new')
            self.assertEqual(p.read_bytes(),b'new');self.assertEqual(list(Path(d).iterdir()),[p])

if __name__=='__main__': unittest.main()
