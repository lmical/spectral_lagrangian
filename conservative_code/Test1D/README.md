# How to run 1D tests

The [source folder](/Src1D) is in the main folder and it must be compiled from the folder [Test1D](/Test1D) (where a Makefile for gfortran is located), so that the bin1D and mod1D folders will be created in Test1D.

In the [makefile](/Test1D/Makefile.1D.gfortran) one can switch between different models commenting and uncommenting the first lines.

Compile as:
```
make -f Makefile.1D.gfortran dec
```


To  run a test go into a Test folder, for example [Test1D/Test1D_Burgers](/Test1D/Test1D_Burgers), edit the don1d file with the desired parameters and then type from that folder
```
../bin1D/main_dec.out
```
or
```
../bin1D/main_dec.out 1000
```
where the number is the number of cells.

The results are saved in the folder TEST inside the Test1D_Burgers folder.

To visualize the solution, one can use the python scripts in folder [Test1D/viewerScripts1D/](/Test1D/viewerScripts1D/)
where scripts are available for systems of 1 [viewer1d_scalar.py](/Test1D/viewerScripts1D/viewer1d_scalar.py), 2 [viewer1dsw.py](/Test1D/viewerScripts1D/viewer1dsw.py) and 3 equations [viewer1d.py](/Test1D/viewerScripts1D/viewer1d.py). Run the script specifying the TEST folder at the end: example in folder [Test1D/Test1D_Burgers](/Test1D/Test1D_Burgers) run

```
python ../viewerScripts1D/viewer1d_scalar.py --dir TEST
```

One can use the arrows to change between different timesteps.








