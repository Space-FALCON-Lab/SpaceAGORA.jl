import contextlib,hashlib,importlib.util,io,json,math,pathlib,sys,tempfile,unittest
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
from PIL import Image
ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.dont_write_bytecode=True
def load(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/'scripts/dev/terrain'/f'{name}.py')
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
f=load('fetch_moon_site');c=load('check_site_registration');stub=load('stub_quadtree_site')

class TerrainPipeline(unittest.TestCase):
    def test_lola_seam_and_bounds(self):
        # An independent global raster, with different row and column gradients.
        arr=(np.arange(18)[:,None]*100+np.arange(36)[None,:]).astype('<i2')
        calls=[]
        def fetch(url,headers):
            a,b=map(int,headers['Range'][6:].split('-'));calls.append((a,b))
            self.assertEqual(a//72,b//72)
            return arr.tobytes()[a:b+1]
        with patch.multiple(f,LOLA_PPD=.1,LOLA_LINES=18,LOLA_SAMPLES=36,LOLA_RECORD_BYTES=72),patch.object(f,'fetch',fetch):
            for left,right in ((-20.,20.),(340.,380.),(-380.,-340.)):
                g,s,n,w,e=f.lola_window(20,-20,left,right,1,1,'test')
                np.testing.assert_array_equal(g,arr[7:11][:,[34,35,0,1]]*.5)
                self.assertEqual((s,n,w,e),(-20,20,left,right))
            g,s,n,w,e=f.lola_window(-71,-90,340,380,2,2,'pole')
            np.testing.assert_array_equal(g,np.array([[arr[16:18,34:36].mean()*.5,arr[16:18,0:2].mean()*.5]]))
            self.assertEqual((s,n),(-90,-70))
            for values in [(float('nan'),0,1,2),(1,0,2,1),(91,0,1,2),(1,0,0,360)]:
                with self.assertRaises(ValueError):f.lola_window(*values,1,1,'invalid')
        self.assertTrue(calls)

    def test_grid_metadata_and_reuse(self):
        with tempfile.TemporaryDirectory() as d:
            base=pathlib.Path(d)/'dem';a=np.array([[1,2,3],[4,5,6]],dtype=float)
            meta=f.write_grid(base,a,-2,2,170,190,'synthetic')
            self.assertEqual((meta['rows'],meta['cols']),(2,3));self.assertTrue(f.existing(base))
            self.assertEqual(meta['reference_radius_m'],1737400)
            np.testing.assert_array_equal(np.fromfile(base.with_suffix('.f32'),dtype='<f4').reshape(2,3),a)
            with self.assertRaises(ValueError):f.existing(base,{'different':'request'})
            base.with_suffix('.f32').write_bytes(b'bad')
            with self.assertRaises(ValueError):f.existing(base)
            for a in (np.array([[np.nan]]),np.empty((0,2)),np.array([1,2])):
                with self.assertRaises(ValueError):f.write_grid(base,a,-2,2,170,190,'bad')

    def test_nac_centers_survive_stride(self):
        # Synthetic real GeoTIFF with area pixels, north-up lunar projection.
        double=(1737400.,1.,180.)
        keys=[1,1,0,8,1024,0,1,1,1025,0,1,1,3075,0,1,17,3076,0,1,9001,2057,34736,1,0,2058,34736,1,0,3078,34736,1,1,3080,34736,1,2]
        scale=(2.,2.,0.);tie=(0.,0.,0.,-6.,6.,0.)
        geo=f.equirectangular_geometry(scale,tie,keys,double)
        with tempfile.TemporaryDirectory() as d:
            out=pathlib.Path(d);(out/'raw').mkdir();a=(np.arange(35).reshape(5,7)-2000).astype(np.float32)
            raw=out/'raw'/'nac.tif';Image.fromarray(a).save(raw,tiffinfo={33550:scale,33922:tie,34735:tuple(keys),34736:double})
            spec={'url':'https://example.invalid/nac.tif','lines':5,'samples':7}
            with patch.dict(f.NAC_DTM,synthetic=spec):
                meta=f.fetch_nac(out,'synthetic',0,180,.001,4,None)
            sub=np.fromfile(out/'dem_nac.f32',dtype='<f4').reshape(meta['rows'],meta['cols'])
            np.testing.assert_array_equal(sub,a[::2,::2])
            for r in range(sub.shape[0]):
                for col in range(sub.shape[1]):
                    lon=meta['lon_min']+(col+.5)*(meta['lon_max']-meta['lon_min'])/sub.shape[1]
                    lat=meta['lat_max']-(r+.5)*(meta['lat_max']-meta['lat_min'])/sub.shape[0]
                    expected=f.projected_lonlat(geo,2*col+.5,2*r+.5)
                    self.assertAlmostEqual(lon,expected[0],places=11);self.assertAlmostEqual(lat,expected[1],places=11)
            self.assertTrue(np.all(sub<0),'Do not add/subtract lunar radius based on sample magnitude.')
        for bad in (float('nan'),-1):
            with self.assertRaises(ValueError):f.equirectangular_geometry((bad,2,0),tie,keys,double)

    def test_registration_and_corridor(self):
        rng=np.random.default_rng(21);a=rng.normal(size=(64,64));b=np.roll(a,(3,-5),(0,1))
        r,col,peak=c.phase_correlate(a,b);self.assertAlmostEqual(r,-3,delta=.1);self.assertAlmostEqual(col,5,delta=.1);self.assertGreater(peak,.5)
        for bad in (np.zeros((5,5)),np.full((5,5),np.nan)):
            with self.assertRaises(ValueError):c.phase_correlate(bad,bad)
        rows,cols=np.indices((5,7));height=2*rows+3*cols
        shade=c.hillshade(height,(2,3),0,90);np.testing.assert_allclose(shade,1/math.sqrt(3),atol=1e-15)
        root=f.quadtree_root([(0,179,0),(0,-179,0)],3)
        self.assertLess(root['lon_max']-root['lon_min'],23)
        self.assertGreater(root['lon_min'],160)

    def test_range_response_fails_closed(self):
        class Response(io.BytesIO):
            status=200;headers={}
        with patch.object(f.urllib.request,'urlopen',return_value=Response(b'full raster')):
            with self.assertRaises(ValueError):f.fetch('https://example.invalid',headers={'Range':'bytes=2-5'})


    def test_seam_track_and_actual_trek_reads(self):
        wrapped=[(0,179,0),(0,-179,0)]; unwrapped=[(0,179,0),(0,181,0)]
        root=f.quadtree_root(wrapped,3)
        self.assertEqual(f.level_coverage(wrapped,root,5,.2), f.level_coverage(unwrapped,root,5,.2))
        self.assertTrue(f.level_coverage(wrapped[:1],root,5,.2))
        self.assertTrue(f.level_coverage(wrapped[1:],root,5,.2))
        seen=[]
        def tile(url,**kwargs):
            x=int(url.rsplit('/',1)[1].split('.')[0]); seen.append(x)
            b=io.BytesIO();Image.new('L',(256,256),70 if x==1 else 200).save(b,format='PNG');return b.getvalue()
        with tempfile.TemporaryDirectory() as d,patch.object(f,'fetch',tile):
            trek=f.TrekTiles([d]);root=dict(lat_min=-90,lat_max=90,lon_min=90,lon_max=270)
            image,key=f.node_image(trek,{},root,0,0,0,0,256)
            a=np.asarray(image)
            self.assertEqual(key,'wac');self.assertTrue(np.all(a[:,:128]==70));self.assertTrue(np.all(a[:,128:]==200))
            self.assertEqual(seen,[1,0]);self.assertIsNone(trek.get('wac',0,0,-1))
            self.assertEqual(len(seen),2)

    def test_texture_matching_wraps(self):
        with tempfile.TemporaryDirectory() as d:
            path=pathlib.Path(d)/'texture.png'
            a=np.array([[10,20,30,40],[10,20,30,40]],dtype=np.uint8);Image.fromarray(a).save(path)
            root=dict(lat_min=-90,lat_max=90,lon_min=90,lon_max=270)
            source=Image.fromarray(np.array([[40,10],[40,10]],dtype=np.uint8))
            image,(gain,offset)=f.match_globe_texture(source,path,root)
            np.testing.assert_array_equal(image,source);self.assertAlmostEqual(gain,1);self.assertAlmostEqual(offset,0)
            legacy=stub.Source(str(path),-90,90,-180,180,'global')
            self.assertTrue(legacy.contains(-90,90,90,270))
            np.testing.assert_array_equal(np.asarray(legacy.crop(-90,90,90,270,2))[:,:,0],np.asarray(source))

    def test_registration_reports_independent_axis_distances(self):
        class Source:
            def window(self,*args):return Image.new('LA',(args[-1],args[-1]),(100,255))
        source=Source();meta=dict(rows=20,cols=20,lat_min=50,lat_max=70,lon_min=-30,lon_max=30)
        with patch.object(c,'phase_correlate',return_value=(2.,3.,.8)),contextlib.redirect_stdout(io.StringIO()):
            result=c.against_dem(source,meta,np.arange(400).reshape(20,20),60,0,8*f.M_PER_DEG,[90],45)
            self.assertAlmostEqual(result[2],-2*math.pi*1737400/180)
            self.assertAlmostEqual(result[3],3*3*math.pi*1737400/180*.5,places=7)
            with patch.object(c,'trek_window',return_value=np.ones((8,8))):
                n,e,_=c.against_trek(source,[],'wac',4,60,0,[8])
                expected=math.pi*1737400/(16*256)
                self.assertAlmostEqual(n,-2*expected);self.assertAlmostEqual(e,3*expected*.5)

    @staticmethod
    def tile_options():
        return SimpleNamespace(root_zoom=3,max_level=1,tile_px=256,lod_factor=2.,fine_lod_factor=2.5,
            detail_std=3.2,detail_max_gain=2.2,quality_ramp_level=12,quality_coarse=44,quality_fine=66,
            approach_azimuth=270.,uprange_km=10.,match_texture=None,retile=False,reuse=None,workers=1)

    def test_imagery_request_reuse_and_transactional_rebuild(self):
        a=self.tile_options();track=[(0,0,0),(0,1,0)];root=f.quadtree_root(track,3)
        pixels=np.tile(np.arange(256,dtype=np.uint8),(256,1))
        def fake_image(*args):return Image.fromarray(pixels),'wac'
        with tempfile.TemporaryDirectory() as d,patch.object(f.TrekTiles,'prefetch'),patch.object(f,'node_image',fake_image),contextlib.redirect_stdout(io.StringIO()):
            out=pathlib.Path(d);meta=f.build_tiles(out,track,root,a)
            self.assertTrue(meta['nodes']);self.assertEqual(meta['producer_version'],2)
            before={str(x.relative_to(out/'imagery')):x.read_bytes() for x in (out/'imagery').rglob('*') if x.is_file()}
            with patch.object(f,'_build_tiles',side_effect=AssertionError('cache should be reused')):
                self.assertEqual(f.build_tiles(out,track,root,a),meta)
            a.detail_std=2.
            with self.assertRaisesRegex(ValueError,'does not match'):f.build_tiles(out,track,root,a)
            a.retile=True
            with patch.object(f,'node_image',return_value=(None,None)):
                with self.assertRaisesRegex(ValueError,'No root imagery'):f.build_tiles(out,track,root,a)
            self.assertEqual(before,{str(x.relative_to(out/'imagery')):x.read_bytes() for x in (out/'imagery').rglob('*') if x.is_file()})
            with patch.object(f,'_build_tiles',side_effect=RuntimeError('source unavailable')):
                with self.assertRaisesRegex(RuntimeError,'source unavailable'):f.build_tiles(out,track,root,a)
            self.assertEqual(before,{str(x.relative_to(out/'imagery')):x.read_bytes() for x in (out/'imagery').rglob('*') if x.is_file()})
            self.assertFalse(list(out.glob('.imagery-build-*')))
            rebuilt=f.build_tiles(out,track,root,a)
            self.assertEqual(rebuilt['request']['options']['detail_std'],2.)
            a.retile=False;tile=out/'imagery'/rebuilt['nodes'][0]['file'];tile.write_bytes(b'corrupt')
            with self.assertRaisesRegex(ValueError,'checksum mismatch'):f.build_tiles(out,track,root,a)

    def test_archive_cache_checks_url_and_geometry(self):
        with tempfile.TemporaryDirectory() as d:
            source=f.ArchiveRaster('a11_pho_r',d,verbose=False)
            source.spec=dict(source.spec);source.spec.pop('bbox');source.spec.pop('label_m_per_px')
            idx=dict(index_version=2,source_url=source.spec['url'],width=2,height=2,tile_w=2,tile_h=2,
                tiles_across=1,tiles_down=1,tile_bytes=4,x0=0.,y0=0.,sx=1.,sy=1.,radius_m=1737400.,
                std_parallel_deg=0.,center_lon_deg=0.,false_easting_m=0.,false_northing_m=0.,tile_offsets=[16])
            (source.dir/'index.json').write_text(json.dumps(idx));self.assertEqual(source.index()['width'],2)
            source._index=None;idx['source_url']='https://example.invalid/other'
            (source.dir/'index.json').write_text(json.dumps(idx))
            with self.assertRaisesRegex(ValueError,'mismatched'):source.index()
            with self.assertRaisesRegex(ValueError,'mismatched'):source.index()
            source._index=None;idx['source_url']=source.spec['url'];idx['radius_m']=1.
            (source.dir/'index.json').write_text(json.dumps(idx))
            with self.assertRaisesRegex(ValueError,'datum'):source.index()

    def test_http_range_acceptance(self):
        class Response(io.BytesIO):
            status=206;headers={'Content-Range':'bytes 2-5/20'}
        for payload,accepted in [(b'abcd',True),(b'abc',False),(b'abcde',False)]:
            with patch.object(f.urllib.request,'urlopen',return_value=Response(payload)):
                if accepted:self.assertEqual(f.fetch('https://example.invalid',headers={'Range':'bytes=2-5'}),payload)
                else:
                    with self.assertRaisesRegex(ValueError,'range response'):f.fetch('https://example.invalid',headers={'Range':'bytes=2-5'})
        response=Response(b'abcd');response.headers={'Content-Range':'bytes 1-4/20'}
        with patch.object(f.urllib.request,'urlopen',return_value=response):
            with self.assertRaisesRegex(ValueError,'exact byte range'):f.fetch('https://example.invalid',headers={'Range':'bytes=2-5'})

    def test_cli_invalid_profile_precedes_any_download(self):
        with tempfile.TemporaryDirectory() as d,patch.object(f,'fetch',side_effect=AssertionError('network forbidden')),contextlib.redirect_stderr(io.StringIO()):
            out=pathlib.Path(d)/'out';profile=pathlib.Path(d)/'profile.json';profile.write_text('[[0,1],[0,0]]')
            with self.assertRaises(SystemExit):f.main(['--site','0.67416','23.47314','--out',str(out),'--profile',str(profile)])
            self.assertFalse(out.exists())
            with self.assertRaises(SystemExit):f.main(['--site','0.67416','23.47314','--out',str(out),'--detail-max-gain','nan'])
            self.assertFalse(out.exists())

    def test_legacy_stub_build_and_protection(self):
        with tempfile.TemporaryDirectory() as d,contextlib.redirect_stdout(io.StringIO()):
            root=pathlib.Path(d);site=root/'site';site.mkdir();out=root/'stub';texture=root/'texture.png'
            Image.new('RGB',(2048,1024),(80,100,120)).save(texture)
            dem=f.write_grid(site/'dem',np.array([[1.,2.],[3.,4.]]),-1,1,-1,1,'synthetic')
            metadata={'site':{'lat_deg':0.,'lon_deg':0.},'dem':[{'name':'dem',**dem}]}
            (site/'site.json').write_text(json.dumps(metadata))
            args=['--site',str(site),'--out',str(out),'--global-texture',str(texture),'--root-deg','8','--uprange-km','1','--tile-px','8','--max-level','1']
            failout=root/'too_coarse'
            with self.assertRaisesRegex(ValueError,'No root imagery'):
                stub.main(['--site',str(site),'--out',str(failout),'--global-texture',str(texture),'--root-deg','.25','--uprange-km','1','--tile-px','16','--max-level','0'])
            self.assertFalse((failout/'site.json').exists())
            self.assertEqual(stub.main(args),0)
            output=json.loads((out/'site.json').read_text());tiles=json.loads((out/output['tiles']).read_text())
            self.assertTrue(tiles['nodes']);self.assertEqual((out/'dem.f32').read_bytes(),(site/'dem.f32').read_bytes())
            for node in tiles['nodes']:self.assertTrue((out/'tiles'/node['file']).is_file())
            with self.assertRaisesRegex(ValueError,'already exist'):stub.main(args)
            self.assertEqual((out/'dem.f32').read_bytes(),(site/'dem.f32').read_bytes())
            with self.assertRaisesRegex(ValueError,'must differ'):stub.main(['--site',str(site),'--out',str(site)])
            with self.assertRaisesRegex(ValueError,'already has a quadtree'):stub.main(['--site',str(out),'--out',str(root/'new')])

    def test_reused_dems_are_self_contained(self):
        with tempfile.TemporaryDirectory() as d,contextlib.redirect_stdout(io.StringIO()):
            parent=pathlib.Path(d);original=parent/'original';original.mkdir();out=parent/'reused';out.mkdir()
            grid=np.array([[1.,2.],[3.,4.]],dtype=np.float32)
            with patch.object(f,'lola_window',return_value=(grid,-1.,1.,-1.,1.)):
                original_meta=f.fetch_lola(original,0.,0.,1.,None)
            with patch.object(f,'fetch',side_effect=AssertionError('network forbidden')):
                meta=f.fetch_lola(out,0.,0.,1.,str(original))
            self.assertEqual(meta,original_meta)
            for ext in ('.json','.f32'):
                path=out/('dem_lola'+ext)
                self.assertFalse(path.is_symlink());self.assertEqual(path.resolve().parent,out.resolve())
                self.assertEqual(path.read_bytes(),(original/path.name).read_bytes())
            path=out/'dem_lola.f32';path.unlink();path.symlink_to(original/path.name)
            self.assertTrue(f.borrow(out,str(original),path.name));self.assertFalse(path.is_symlink())

    def test_lola_only_cli_produces_complete_bundle_offline(self):
        raster=np.arange(18*36,dtype='<i2').reshape(18,36)
        def fetch(url,headers):
            lo,hi=map(int,headers['Range'][6:].split('-'))
            return raster.tobytes()[lo:hi+1]
        with tempfile.TemporaryDirectory() as d,patch.multiple(f,LOLA_PPD=.1,LOLA_LINES=18,LOLA_SAMPLES=36,LOLA_RECORD_BYTES=72),patch.object(f,'fetch',fetch),contextlib.redirect_stdout(io.StringIO()):
            out=pathlib.Path(d)/'site'
            f.main(['--site','0.67416','23.47314','--out',str(out),'--no-nac','--no-imagery','--uprange-km','1','--lola-wide-samples','2','--lola-half-deg','.2'])
            meta=json.loads((out/'site.json').read_text());self.assertIsNone(meta['tiles'])
            self.assertEqual([item['name'] for item in meta['dem']],['dem_lola','dem_lola_wide'])
            for item in meta['dem']:
                self.assertTrue(f.existing(out/item['name']))
                self.assertEqual(item['reference_radius_m'],1737400.)

if __name__=='__main__':unittest.main(verbosity=2)
