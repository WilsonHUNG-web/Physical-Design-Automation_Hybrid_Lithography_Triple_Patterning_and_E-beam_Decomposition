# Physical-Design-Automation_Hybrid_Lithography_Triple_Patterning_and_E-beam_Decomposition

--how to compile 
In ./src directory:

(1)Makefile:
$ make

(2)Run the cases:
$ ../bin/main [dmin] [../testcase/inputfile] [../output/outputfile]

e.g.
$ ../bin/main 0.0128 ../testcase/M0.out ../output/M0_00128.out

(3)Verify the solution
$ ./Verify  [dmin] [../testcase/inputfile] [../output/outputfile]

e.g. 
$ ./Verify 0.0128 ../testcase/M0.out ../output/M0_00128.out


Ps. 
You can also use the bash shell... runall.sh by the command:
$ bash runall.sh

In the runall.sh, it looks like:

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
