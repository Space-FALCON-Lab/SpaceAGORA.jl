"""Exercise the reusable workflow shell with temporary scripts and a fake Julia command."""
from pathlib import Path
import os, subprocess, tempfile
w=(Path(__file__).resolve().parents[2]/'.github/workflows/contributor-ci.yml').read_text()
def block(name):
 s=w.split('      - name: '+name+'\n',1)[1].split('        run: |\n',1)[1].split('\n      - name:',1)[0]
 lines=[]
 for line in s.splitlines():
  if line and not line.startswith('          '):break
  lines.append(line[10:])
 return '\n'.join(lines)+'\n'
env=os.environ|dict(SPACEAGORA_REV='a'*40,KIND='companion',MANIFEST='core',DEPWARN='error',JULIA_NUM_THREADS='4',LOAD_ORDERS='true',STARTER='false',CORE_FROM_CALLER='false',REPOSITORY='',REPOSITORY_REF='',DIRECTORY='.',WORKFLOW_SHA='a'*40,WORKFLOW_REF='fixture')
checks=0
for value,ok in [('test/runtests.jl',True),('test/contributor_ci.jl',True),('',False),('/tmp/no.jl',False),('../no.jl',False),('test/../../no.jl',False)]:
 r=subprocess.run(['bash','-c',block("Check the inputs and this workflow's commit")],env=env|{'TEST_ENTRYPOINT':value},capture_output=True,text=True)
 assert (r.returncode==0)==ok,(value,r.stderr,r.stdout)
 checks+=1
with tempfile.TemporaryDirectory(prefix='spaceagora-entrypoint-') as td:
 td=Path(td).resolve();root=td/'repository';(root/'test').mkdir(parents=True);stub=td/'bin';stub.mkdir();out=td/'calls'
 (stub/'julia').write_text('#!/bin/bash\nprintf "%s\\n" "$PWD" "${@: -1}" > "$PROBE_OUTPUT"\n');(stub/'julia').chmod(0o755)
 for name in ('runtests.jl','contributor_ci.jl','with space.jl'):(root/'test'/name).write_text('# fixture\n')
 outside=td/'outside.jl';outside.write_text('# outside\n');(root/'test'/'escape.jl').symlink_to(outside)
 for value,ok in [('test/runtests.jl',True),('test/contributor_ci.jl',True),('test/with space.jl',True),('test/missing.jl',False),('test/escape.jl',False)]:
  r=subprocess.run(['bash','-c',block('Tests')],env=env|dict(TEST_ENTRYPOINT=value,CONTRIBUTOR_REPOSITORY=str(root),CONTRIBUTOR_ENVIRONMENT=str(td/'env'),SPACEAGORA_CONFORMANCE_REPORT_DIR=str(td/'reports'),PATH=str(stub)+':'+os.environ['PATH'],PROBE_OUTPUT=str(out)),capture_output=True,text=True)
  assert (r.returncode==0)==ok,(value,r.stderr,r.stdout)
  if ok:assert out.read_text().splitlines()==[str(root/'test'),str(root/value)]
  checks+=1
print(f'{checks} actual workflow shell cases passed; default/alternate/quoted paths, cwd, traversal and symlink refusals.')
