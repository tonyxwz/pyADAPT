#!/usr/bin/env python
"""
find out which time step in which iteration caused the stuck
"""
import numpy as np
import re
import sys


def whosefault(log_path):

    with open(log_path, "r") as log_file:
        logstr = log_file.read()

    reg = re.compile(r"iter:\s*(\d*)\s*,\s*ts:\s*(\d*)")
    res = reg.findall(logstr)

    reg_nts = re.compile(r"n_ts:\s*(?P<n_ts>\d*)")
    n_ts = reg_nts.search(logstr)
    if n_ts is not None:
        n_ts = int(n_ts["n_ts"])

    reg_niter = re.compile(r"n_iter:\s*(?P<n_iter>\d*)")
    n_iter = reg_niter.search(logstr)
    if n_iter is not None:
        n_iter = int(n_iter["n_iter"])

    counters = np.zeros((n_iter, n_ts), dtype=np.int)

    for match in res:
        iter = int(match[0])
        ts = int(match[1])
        counters[iter, ts] += 1

    reg_init = re.compile(r"init iter (?P<iter>\d+)")
    inits = reg_init.findall(logstr)
    for match in inits:
        iter = int(match)
        counters[iter, 0] += 1

    print("iter\tts\tpid\tattempt\twhen")
    for iter in range(len(counters)):
        if counters[iter, 0] != 0 and counters[iter, -1] == 0:
            row = counters[iter, 1:]
            # not reached the last time step yet
            # safely assume that steps before last 1 flag are finished
            pos = (row != 0).argmin()

            # fixme: if an iter misses one record and not finished
            # if np.all(row[pos[-1]+1:] == 0):
            when = re.search(
                r"(?P<date>\d*-\d*-\d*)\s*(?P<time>\d*:\d*:\d*,\d*).*"
                + r"\s+(?P<pid>\d+)\s+"
                + r"iter:\s*"
                + str(iter)
                + r"\s*,\s*ts:\s*"
                + str(pos),
                logstr,
            )
            # when = reg2.search(logstr)
            if when is None:
                # if iter is not run yet
                continue
            # iteration () dies at timestep ()
            # FIXME: row[0] is containing attempt_limit+1
            print(
                f"{iter}\t{pos}\t{when['pid']}\t{row[0]}\t{when['date']} {when['time']}"
            )


def main():
    if len(sys.argv) > 1:
        log = sys.argv[1]
    else:
        import glob

        logs = glob.glob("adapt_" + "*.log")

        latest = 0.0
        for _log in logs:
            date, time = _log.rsplit(".", 1)[0].rsplit("_")[-2:]
            time = float("".join(date.split("-") + time.split(".")))
            if time > latest:
                latest = time
                log = _log
    print(log)
    whosefault(log)


if __name__ == "__main__":
    # whosefault("adapt_2020-05-20_02.06.04.log", n_iter=120, n_ts=169)
    main()
