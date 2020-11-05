from time import time
from pyKVFinder.modules.utils import read_vdw_dat, read_vdw_json, read_vdw_toml

def run():

    print("\n[===> Reading van der Waals radii")

    print("\n> from vdw.dat")
    start = time()
    vdw = read_vdw_dat("data/vdw.dat")
    end = time()
    dat_time = end-start
    print(vdw)
    
    print("\n> from vdw.toml")
    start = time()
    vdw = read_vdw_toml("data/vdw.toml")
    end = time()
    toml_time = end-start
    print(vdw)

    print("\n> from vdw.json")
    start = time()
    vdw = read_vdw_json("data/vdw.json")
    end = time()
    json_time = end-start
    print(vdw)

    print("\nLoading Time")
    print(f"DAT:\t {dat_time:0.6f}")
    print(f"TOML:\t {toml_time:0.6f}")
    print(f"JSON:\t {json_time:0.6f}")