"""Run one serial experiment process with recorded background-process checks."""
import datetime,json,os,signal,subprocess,sys,time
from pathlib import Path
path=Path(sys.argv[1]);command=sys.argv[2:]
if path.exists():raise SystemExit("Refusing to overwrite process monitor evidence")
if not command:raise SystemExit("usage: run-performance-process.py NEW_MONITOR_JSON COMMAND...")
def snapshot():
    raw=subprocess.check_output(["ps","-axo","pid=,ppid=,state=,time=,rss=,command="],text=True)
    rows=[];parents={}
    for line in raw.splitlines():
        fields=line.strip().split(None,5)
        if len(fields)!=6:continue
        pid,ppid,state,cpu,rss,args=fields
        parents[int(pid)]=int(ppid)
        if os.path.basename(args.split()[0])=="julia":
            rows.append(dict(pid=int(pid),ppid=int(ppid),state=state,cpu_time=cpu,rss_kib=int(rss),command=args))
    return dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),julia=rows),parents
def descendant(pid,parent,parents):
    seen=set()
    while pid and pid not in seen:
        if pid==parent:return True
        seen.add(pid);pid=parents.get(pid,0)
    return False
initial,_=snapshot();allowed={}
for item in initial["julia"]:
    if item["pid"]!=55147 or not item["state"].startswith("S"):
        raise SystemExit("Another active/unknown Julia process is present: "+str(item))
    allowed[item["pid"]]=item["cpu_time"]
memory={}
for key,cmd in [("vm_stat",["vm_stat"]),("swap",["sysctl","vm.swapusage"]),("pressure",["memory_pressure","-Q"])]:
    try:memory[key]=subprocess.check_output(cmd,text=True,stderr=subprocess.STDOUT,timeout=10)
    except Exception as exc:memory[key]=str(exc)
record=dict(status="before execution",command=command,initial=initial,documented_idle_cpu=allowed,memory_context=memory,observations=[],monitor_interval_s=2.)
def persist():
    path.write_text(json.dumps(record,indent=2)+"\n")
    compact=["status = "+json.dumps(record["status"]),"exit_code = "+str(record.get("exit_code",-999)),"command = "+json.dumps(command),"observations = "+str(len(record["observations"])),"background_activity_valid = "+str("violation" not in record).lower()]
    path.with_suffix(".toml").write_text("\n".join(compact)+"\n")
persist();process=subprocess.Popen(command);record["pid"]=process.pid;violation=None
while True:
    observation,parents=snapshot();record["observations"].append(observation)
    for item in observation["julia"]:
        if descendant(item["pid"],process.pid,parents):continue
        if item["pid"] not in allowed or item["cpu_time"]!=allowed[item["pid"]] or not item["state"].startswith("S"):
            violation="Non-idle/unknown Julia observed: "+str(item)
    if violation:
        record["violation"]=violation
        if process.poll() is None:process.send_signal(signal.SIGINT)
        process.wait();record["status"]="invalidated by background activity";record["exit_code"]=process.returncode;persist();raise SystemExit(3)
    code=process.poll()
    if code is not None:
        record["status"]="terminal clean observations" if code==0 else "child failed; no performance verdict"
        record["exit_code"]=code;persist();raise SystemExit(code if code>=0 else 128-code)
    persist();time.sleep(2)
