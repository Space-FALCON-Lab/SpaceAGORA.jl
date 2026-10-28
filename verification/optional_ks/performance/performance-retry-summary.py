"""Additive postprocessing; original frozen experiment programs stay unchanged.
Requires Python3.11+. Run after all selected workers are terminal, before changing
installed sources. Explicit selection is by valid chronology, never measured speed.
"""
import datetime,hashlib,json,math,statistics,sys,tomllib
from pathlib import Path

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def toml(path):return tomllib.loads(path.read_text())
def check(condition,message):
    if not condition:raise RuntimeError(message)
def finite_positive(value):return isinstance(value,(int,float)) and math.isfinite(value) and value>0

def summarize(manifest_path,out):
    root=manifest_path.resolve().parent
    manifest=json.loads(manifest_path.read_text())
    check(not out.exists(),'Refusing to overwrite summary')
    def artifact(entry,key):
        path=root/entry[key]
        check(digest(path)==entry[key+'_sha256'],'Artifact hash mismatch: '+str(path))
        return path
    plan=json.loads(artifact(manifest,'selection_plan').read_text())
    expected=[(c,m,str(Path(d)/(c+'-'+m))) for c,m,d in plan['workers']]
    actual=[(e['case'],e['method'],str(Path(e['timing']).parent)) for e in manifest['workers']]
    check(actual==expected,'Selection differs from prospective plan')
    check([(e['reason'],[f['path'] for f in e['files']]) for e in manifest['excluded']]==[(e['reason'],e['paths']) for e in plan['excluded']],'Exclusion ledger differs from plan')
    cal=toml(artifact(manifest,'calibration'))
    check(cal['status']=='matched accuracy selected; no speed measurements','Calibration incomplete')
    identity=cal['identity'];cal_sha=manifest['calibration_sha256']
    for name,sha in identity['scripts'].items():check(digest(root/name)==sha,'Original experiment script drift: '+name)
    results={};seen=set();hardware=None
    for entry in manifest['workers']:
        case,mode=entry['case'],entry['method'];pair=(case,mode)
        check(case in ('two_body','J2_plus_drag') and mode in ('ks','dp8') and pair not in seen,'Unexpected/duplicate selection')
        seen.add(pair)
        r=toml(artifact(entry,'timing'));mt=toml(artifact(entry,'monitor_toml'));mj=json.loads(artifact(entry,'monitor_json').read_text())
        check(mt['status']==mj['status']=='terminal clean observations' and mt['exit_code']==mj['exit_code']==0,'Worker not terminal clean0')
        check(mt['background_activity_valid'] and 'violation' not in mj,'Invalidated worker')
        check(mt['command']==mj['command'] and len(mj['observations'])==mt['observations']>0,'Monitor mismatch')
        cmd=mt['command'];launch=Path(entry['launch_working_directory'])
        check(cmd[-3]==case and cmd[-2]==mode,'Misassociated command case/method')
        resolved=lambda p:(launch/Path(p)).resolve()
        check(resolved(cmd[-1])==resolved(str(Path(entry['timing']).parent)),'Misassociated command output path')
        check(resolved(cmd[-4])==resolved(manifest['calibration']),'Misassociated command calibration path')
        timestamps=[datetime.datetime.fromisoformat(o['utc']) for o in mj['observations']]
        gaps=[(b-a).total_seconds() for a,b in zip(timestamps,timestamps[1:])]
        coverage={'maximum_observation_gap_s':max(gaps,default=0.),'gaps_over_10s':sum(g>10 for g in gaps),'note':'Post-experiment timing-quality diagnostic, not a retroactive change to numerical acceptance or the original process policy. Long unobserved intervals preclude a reliable wall-time ratio.'}
        allowed=mj['documented_idle_cpu'];own=mj['pid']
        for obs in mj['observations']:
            parents={x['pid']:x['ppid'] for x in obs['julia']}
            for p in obs['julia']:
                if str(p['pid']) in allowed:
                    check(p['state'].startswith('S') and p['cpu_time']==allowed[str(p['pid'])],'Background counter/state changed')
                else:
                    pid=p['pid'];visited=set()
                    while pid not in visited and pid!=own:
                        visited.add(pid);pid=parents.get(pid,0)
                    check(pid==own,'Unknown Julia process in selected observations')
        check(r['case']==case and r['method']==mode and r['identity']==identity and r['calibration_sha256']==cal_sha,'Worker source/calibration mismatch')
        check(r['status']=='complete worker; paired accuracy summary still required','Incomplete worker')
        check(hardware is None or hardware==r['hardware'],'Mixed hardware');hardware=r['hardware']
        settings={'rtol':1e-12,'atol':1e-13,'dtmax_s':5.} if mode=='ks' else {k:cal['selected'][case][k] for k in ('rtol','atol','dtmax_s')}
        check(all(r[k]==v for k,v in settings.items()),'Worker settings differ from calibrated selection')
        target=cal['KS_current'][case]['accuracy'];runs=r['runs']
        check(runs[0]['kind']=='first use' and runs[0]['index']==0 and all(x['kind']=='warm candidate' and x['index']==i for i,x in enumerate(runs[1:],1)),'Run order mismatch')
        check(3<=len(runs)-1<=5,'Warm call count outside policy')
        stable=[]
        for row in runs:
            c=row['contract'];a=row['accuracy']
            check(c['complete'] and c['endpoint_s']==86400. and c['retcode']=='Success','Physical endpoint/retcode gate')
            check(c['clock_applicable']==(mode=='ks'),'Clock applicability mismatch')
            if mode=='ks':check(math.isfinite(c['clock_residual_s']) and 0<=c['clock_residual_s']<=1e-7,'Physical clock gate')
            check(all(math.isfinite(v) and v>=0 for v in a.values()),'Invalid raw accuracy')
            common=a['position_m']<=1 and a['velocity_m_s']<=.001 and a['mass_kg']<=1e-9 and a['heat_budget_ratio']<=1
            match=common and all(finite_positive(a[k]) and finite_positive(target[k]) and .5<=a[k]/target[k]<=2 for k in ('position_m','velocity_m_s'))
            check(match and row['matched_accuracy'],'Raw common/calibrated accuracy gate')
            check(finite_positive(row['wall_s']) and all(math.isfinite(row[k]) and row[k]>=0 for k in ('compile_s','recompile_s','gc_s','allocated_bytes')),'Invalid timing counters')
            expected_stable=row['kind']=='warm candidate' and row['compile_s']<=.01*row['wall_s']
            check(row.get('stable_execution',False)==expected_stable,'Residual compilation classification mismatch')
            if expected_stable:stable.append(row)
        check(len(stable)==3,'Exactly3 stable calls required')
        wall=[x['wall_s'] for x in stable]
        results.setdefault(case,{'methods':{}})['methods'][mode]={'monitor_coverage':coverage,'settings':settings,'median_wall_s':statistics.median(wall),'min_wall_s':min(wall),'max_wall_s':max(wall),'stable_runs':stable,'first_use_and_warmups':[x for x in runs if x not in stable],'package_load':r['package_load'],'hardware':r['hardware']}
    check(seen=={(c,m) for c in ('two_body','J2_plus_drag') for m in ('ks','dp8')},'Four selected workers required')
    for case,item in results.items():
        methods=item['methods'];ratios={}
        for key in ('position_m','velocity_m_s'):
            a=[x['accuracy'][key] for x in methods['ks']['stable_runs']];b=[x['accuracy'][key] for x in methods['dp8']['stable_runs']]
            worst=max(max(a)/min(b),max(b)/min(a));check(math.isfinite(worst) and worst<=2,'Actual paired error differs by more than2x')
            ratios[key]=worst
        item['worst_cross_pair_error_ratio']=ratios
        continuous=all(m['monitor_coverage']['gaps_over_10s']==0 for m in methods.values())
        item['performance_verdict']='bounded observed-machine timing; OS/memory/repetition limitations apply' if continuous else 'inconclusive wall-time comparison: major gaps in monitor observations'
        item['KS_over_DP8_median_wall_ratio']=methods['ks']['median_wall_s']/methods['dp8']['median_wall_s'] if continuous else None
    for exclusion in manifest['excluded']:
        check(bool(exclusion['reason']),'Missing exclusion reason')
        for file in exclusion['files']:artifact(file,'path')
    report={'status':'complete bounded matched-accuracy comparison','identity':identity,'calibration_sha256':cal_sha,'cases':results,'excluded':manifest['excluded'],'postprocessing':{'source_sha256':digest(Path(__file__)),'selection_manifest_sha256':digest(manifest_path),'selection_plan_sha256':manifest['selection_plan_sha256'],'note':'Additive retry-aware analysis; not retroactively included in original seven-script experiment identity. Strict original gates retained and raw gates rechecked.'}}
    out.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('ACTUAL_PAIRED_ACCURACY_RETRY_SUMMARY_SUCCESS')
if __name__=='__main__':
    check(len(sys.argv)==3,'usage: performance-retry-summary.py SELECTION.json NEW_SUMMARY.json')
    summarize(Path(sys.argv[1]),Path(sys.argv[2]))
