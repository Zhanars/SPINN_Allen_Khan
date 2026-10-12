"""Wait until a log file contains a marker, then run a job group: python wait_and_run.py <log> <marker> <group> <nproc>"""
import sys, time, subprocess, os
log, marker, group, nproc = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
while True:
    if os.path.exists(log) and marker in open(log, encoding="utf-8", errors="ignore").read():
        break
    time.sleep(30)
subprocess.run([sys.executable, "r1_run.py", group, nproc], stdout=open(os.path.join("..", "data_r1", group + ".log"), "w"), stderr=subprocess.STDOUT)
