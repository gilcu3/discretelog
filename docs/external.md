# External tools

For big numbers a pure python implementation has no chance to finish in a 
reasonable time. Therefore, for factoring really large numbers this
tool uses `yafu` and `cado-nfs` instead of `primefac`, and for computing
discrete logs in fields larger than `10^27` we use `cado-nfs`.

## Yafu

[yafu](https://github.com/bbuhrow/yafu/) is a state of the art hobbist tool to 
factor integers. `discretelog` automatically uses it when the number to factor 
is beyond what `primefac` is able to handle efficiently. 

I personally have the following script on my `PATH`, making `yafu` accesible 
to other programs. I compiled `yafu` on the `~/opt/yafu/` folder, and use it 
with 8 threads. Given that currently the tool generates a lot of temporary 
files, I launch it always from a fresh folder under `/tmp`.

```bash
# ~/bin/yafu
tmp=$(mktemp -d --tmpdir "yafu-XXXXXXXXXX")
echo "Using tmp dir $tmp"
cd $tmp
cp ~/opt/yafu/yafu.ini .
~/opt/yafu/yafu -lathreads 8 -threads 8 $@
```

## Cado-nfs
[cado-nfs](https://gitlab.inria.fr/cado-nfs/cado-nfs) is a state of the art 
academic tool capable of factoring integers and computing the discrete 
logarithms in prime fields. It does not work with small values nor non-prime
moduli, therefore requires something like `discretelog` to be usable in all
cases.

I personally have the following script on my `PATH`, making `cado-nfs` 
accesible to other programs. `cado-nfs` was compiled in `~/opt/cado-nfs`, and
uses all cores available by default.

```bash
# ~/bin/cado-nfs
cd ~/opt/cado-nfs
python cado-nfs.py $@
```
