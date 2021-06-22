# Physical-Design-Automation_Hybrid_Lithography_Triple_Patterning_and_E-beam_Decomposition

--how to compile <br>
In ./src directory:<br>

(1)Makefile:<br>
```
$ make
```
(2)Run the cases:<br>
```
$ ../bin/main [dmin] [../testcase/inputfile] [../output/outputfile]<br>
```
e.g.<br>
```
$ ../bin/main 0.0128 ../testcase/M0.out ../output/M0_00128.out
```
(3)Verify the solution
```
$ ./Verify  [dmin] [../testcase/inputfile] [../output/outputfile]
```
e.g. <br>
```
$ ./Verify 0.0128 ../testcase/M0.out ../output/M0_00128.out
```

Ps. <br>
You can also use the bash shell... runall.sh by the command:<br>
```
$ bash runall.sh
```
In the runall.sh, it looks like:<br>
```
echo "[dmin = 0.0128]"
../bin/main 0.0128 ../testcase/M0.out ../output/M0_0.0128.out
../bin/main 0.0128 ../testcase/aes.out ../output/aes_0.0128.out
echo "     [M0, 0.0128]"
./Verify 0.0128 ../testcase/M0.out ../output/M0_0.0128.out
echo "     [aes, 0.0128]"
./Verify 0.0128 ../testcase/aes.out ../output/aes_0.0128.out
echo "[dmin = 0.064]"
../bin/main 0.064 ../testcase/M0.out ../output/M0_0.064.out
../bin/main 0.064 ../testcase/aes.out ../output/aes_0.064.out
echo "     [M0, 0.064]"
./Verify 0.064 ../testcase/M0.out ../output/M0_0.064.out
echo "     [aes, 0.064]"
./Verify 0.064 ../testcase/aes.out ../output/aes_0.064.out
```
Referenced paper:
Haitong Tian, Hongbo Zhang, Zigang Xiao, Martin D. F. Wong, "Hybrid lithography for triple patterning decomposition and E-beam lithography," Proc. SPIE 9052, Optical Microlithography XXVII, 90520P (31 March 2014); https://doi.org/10.1117/12.2046499
