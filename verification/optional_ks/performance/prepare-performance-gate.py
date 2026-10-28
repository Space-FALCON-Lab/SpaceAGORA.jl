from pathlib import Path
import json,hashlib,subprocess,sys
if len(sys.argv)!=4:raise SystemExit('usage: prepare-performance-gate.py EVIDENCE_ROOT SPACEAGORA_ROOT KS_INSTALLED_ROOT')
D,SG,KR=map(lambda s:Path(s).resolve(),sys.argv[1:])
old=json.loads((D/'adapter-candidate-qualification-source.json').read_text())
for f,sha in old['spaceagora_files_sha256'].items():
 if f.startswith(('src/','ext/')) or f=='Project.toml':
  assert hashlib.sha256((SG/f).read_bytes()).hexdigest()==sha,f
installed=json.loads((D/'complete-source-installed-identity.json').read_text())
qualified=json.loads((D/'ks-complete-candidate-test-source-identity.json').read_text())
assert installed['candidate']==qualified['complete_candidate']=='05d969c24a6595c2752f7ff947af77c504053623'
for f,sha in qualified['source_sha256'].items():assert hashlib.sha256((KR/f).read_bytes()).hexdigest()==sha,f
sources={}
for name,root in [('SpaceAGORA',SG),('KSPropagator',KR)]:
 for base in ['src']+(['ext'] if name=='SpaceAGORA' else []):
  for f in sorted((root/base).rglob('*')):
   if f.is_file():sources[name+'/'+str(f.relative_to(root))]=hashlib.sha256(f.read_bytes()).hexdigest()
 sources[name+'/Project.toml']=hashlib.sha256((root/'Project.toml').read_bytes()).hexdigest()
evidence=['baseline-attempt1/results.toml','adapter-candidate-qualification-source.json','ks-complete-candidate-test-source-identity.json','complete-source-installed-identity.json','logs/complete-source-installed-checks-attempt1.log','review/adapter-current-math-review.json','review/error-display-io-binding-review.json','review/current-upstream-lifecycle-independent-review-attempt1.md','logs/github-checkout-without-ks-checks-attempt1.log']
lines=['status = "correctness prerequisites verified; performance not yet measured"','spaceagora_scientific_commit = "c091191a1ae0dc602ae50025889d34d4b2800a96"','ks_installed_commit = "05d969c24a6595c2752f7ff947af77c504053623"','scientific_source_note = "Original trajectory records retain KS8cfc936 provenance; independently approved Base.IO annotation-only delta and complete4170 suite plus installed316 checks qualify current numerical source. Benchmark calibration also recomputes KS fine errors on current installed source."','[physical_source_sha256]']
lines += [json.dumps(k)+' = '+json.dumps(v) for k,v in sorted(sources.items())]
lines += ['[evidence_sha256]']+[json.dumps(k)+' = '+json.dumps(hashlib.sha256((D/k).read_bytes()).hexdigest()) for k in evidence]
# Run only before calibration. Never change this gate during an experiment.
(D/'performance-correctness-gate.toml').write_text('\n'.join(lines)+'\n')
