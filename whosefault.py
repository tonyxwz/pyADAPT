"""
find out which time step in which iteration caused the stuck
"""
import numpy as np
import re
import sys

def whosefault(log_path, n_iter=None, n_ts=None):

    with open(log_path, "r") as log_file:
        logstr = log_file.read()

    reg = re.compile(r"iter:\s*(\d*)\s*,\s*ts:\s*(\d*)")
    res = reg.findall(logstr)

    flags = np.zeros((n_iter, n_ts))
    for match in res:
        iter = int(match[0])
        ts = int(match[1])
        flags[iter, ts] = 1
    print("iter\tts\twhen")
    for iter in range(len(flags)):
        row = flags[iter, 1:]
        if row[-1] == 0:  # not reached the last time step yet
            # safely assume that steps before last 1 flag are finished
            pos = []
            pos.append((row == 1).argmin())
            # fixme: if an iter misses one record and not finished
            # if np.all(row[pos[-1]+1:] == 0):
            reg2 = re.compile(r"(?P<date>\d*-\d*-\d*)\s*(?P<time>\d*:\d*:\d*,\d*).*iter:\s*"
                                +str(iter)+r"\s*,\s*ts:\s*"
                                +str(pos[-1]))
            when = reg2.search(logstr)
            # iteration () dies at timestep ()
            print(f"{iter}\t{pos}\t{when['date']} {when['time']}")

def main():
    whosefault(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))

if __name__ == "__main__":
    # whosefault("adapt_2020-05-20_02.06.04.log", n_iter=120, n_ts=169)
    main()