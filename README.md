# cpm
A Linux version of CPM flow

Dependencies:

```
sudo apt-get install liblapacke-dev
sudo apt-get install libopencv-dev
```

First, use cpm to get semi dense matches:

```
./cpm img1Name img2Name outMatchName <c> <d> <r> <k> <n>
```

```
c: forward-backward check flag <default: 1>
d: grid spacing <default: 3>
r: search radius <default: 4>
k: pyramid levels <default: 5>
n: iteration times <default: 6>
```

Then, refine the matches to flow:

```
./Match2Flow matches.txt <w> <h>
```

```
w: width of image <no default>
h: height of image <no default>
```