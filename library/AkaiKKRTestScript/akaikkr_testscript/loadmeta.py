# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

import sys
import os
import json
import datetime
import platform
import subprocess

def make_meta(compilermeta = "../meta.json", commitfile="../_commit_id_", cpuinfofile="/proc/cpuinfo"):
    """make minimum meta information automatically."""
    meta = {}
    if os.path.isfile(compilermeta):
        with open(compilermeta) as f:
            _meta = json.load(f)
        meta.update(_meta) 

    if os.path.isfile(commitfile):
        with open(commitfile) as f:
            data = f.read().splitlines()
        s = data[0].split(" ")
        commit = s[1]
        meta["commit"] = commit

    now = datetime.datetime.utcnow().ctime()
    meta["date"] = now

    user = os.environ.get("USER")
    meta["user"] = user

    hostname = os.uname()[1]
    meta["hostname"] = hostname

    _platform = platform.platform()
    meta["platform"] = _platform

    if os.path.isfile(cpuinfofile):
        cmd = ["cat", cpuinfofile]
        try:
            out = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        except:
            pass
        lines =  out.communicate()
        lines = lines[0].decode().splitlines()
        processorname = None
        for x in lines:
            if x.startswith("model name"):
                s = x.split("model name")[1]
                t = s.split(":")[1].strip()
                processorname= t.strip()
                out.stdout.close()
                break
        meta["processorname"] = processorname
    else:
        processor = platform.processor()
        meta["processor"] = processor

    return meta

if __name__ == "__main__":
    print(make_meta())

