# Las unidades son en el SI

# Iniciar
wipe
wipeAnalysis
file mkdir Results;

# -------------------------------------------
#			Generación del modelo
# -------------------------------------------

model basic -ndm 2 -ndf 3 ;#BasicBuilder

# Load of parameters
source Data/Data_HDRB.tcl

# Element height
set h 34.3

# Creation of the element
node 1  0.0 0.0
node 2  0.0 $h

# Restraints
fix  1   1 1 1 
fix  2   0 0 1 

# Axial properties
set fas 3.0
set Kb 461998.9196
set Kr 461998.9196
set Do 70.0
set Di 10.0

# Parameters for cavitation
set Tr 20.40
set a 1.0
set k  20.0
set pm2 0.75
set fc 56.55
set Kvi 2438.284

set bet  0.90;
set eta  0.32;

# Properties of the element
uniaxialMaterial Elastic 5 $Kvi
uniaxialMaterial Elastic 6 $Kb
uniaxialMaterial HDR_mat 2 $a1 $a2 $a3 $fs1 $ps1 $fs2 $ps2 $fs3 $ps3 $fm $pm $ky $fy $bet $eta $pmx $Pfi $Tr
element twoNodeLink 1 1 2 -mat 5 2 6 -dir 1 2 3

puts "Model ok"

timeSeries Constant 7 -factor 1
set pa 345.0
pattern Plain 1 7 {
	load 2 0.0 [expr -$pa] 0.0
}
puts "load ok"

constraints Transformation
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-15 20;# tolerancia #iteraciones
algorithm Newton
integrator LoadControl 1.0;
analysis Static

puts "Before Gravity analysis"

set ok [analyze 1]

puts "Gravity analysis ok"

# Set the gravity loads to be constant & reset the time in the domain
loadConst -time 0.0


# //////////////////////////////////////////////
#          Pseudo-static analysis
# //////////////////////////////////////////////

set IDctrlNode 2
set IDctrlDOF 1
set Tol 1.0e-7;
set LunitTXT "cm"
set Dincr 10.0



# Set lateral load pattern with a Linear TimeSeries
set Plateral 1.0;		# Reference lateral load	
pattern Plain 2 "Linear" {

	load $IDctrlNode $Plateral 0.0 0.0
}

set load_step 1;

# Información a guardar
recorder Element -file ./Results/Displ_HDRB_D1.txt -time -closeOnWrite -ele 1 localDisplacement
recorder Element -file ./Results/Force_HDRB_D1.txt -time -closeOnWrite -ele 1 localForce

set clockI [clock milliseconds];

# set up analysis parameters
source Data/LibAnalysisStaticParameters.tcl;	# constraintsHandler,DOFnumberer,system-ofequations,convergenceTest,solutionAlgorithm,integrator

#  perform Static Cyclic Displacements Analysis
source Data/Displ_D1.tcl

set fmt1 "%s Cyclic analysis: CtrlNode %.3i, dof %.1i, Disp=%.4f %s";	# format for screen/file output of DONE/PROBLEM analysis

set zeroD 0
set D0 0.0 
foreach Dstep $iDstep {
	set D1 $Dstep
	set Dincr [expr ($D1 - $D0)]
	integrator DisplacementControl  $IDctrlNode $IDctrlDOF $Dincr
	analysis Static
	algorithm ModifiedNewton
	# ----------------------------------------------first analyze command------------------------
	set ok [analyze 1]
	# ----------------------------------------------if convergence failure-------------------------
	if {$ok != 0} {
		# if analysis fails, we try some other stuff
		# performance is slower inside this loop	global maxNumIterStatic;# max no. of iterations performed before "failure to converge" is ret'd
		if {$ok != 0} {
			puts "Trying Newton with Initial Tangent .."
			test NormDispIncr   $Tol 1000 0
			algorithm Newton -initial
			analysis Static
			set ok [analyze 1]
			test $testTypeStatic $TolStatic      $maxNumIterStatic    0
			algorithm $algorithmTypeStatic 500 1
		}
		if {$ok != 0} {
			puts "Trying Broyden .."
			algorithm Broyden 24
			analysis Static
			set ok [analyze 1 ]
			algorithm $algorithmTypeStatic
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			algorithm NewtonLineSearch .8
			analysis Static
			set ok [analyze 1 ]
			algorithm $algorithmTypeStatic
		}
		if {$ok != 0} {
			set putout [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
			puts $putout
			return -1
		}; # end if
	}; # end if
	# -----------------------------------------------------------------------------------------------------
	set D0 $D1;			# move to next step
	
	# print load step on the screen
	puts "Load Step: [expr $load_step]"
	set res [expr $load_step%1000]
	
	if {$res == 0} {
		puts "Load Step: [expr $load_step]"
	};
	set load_step [expr $load_step+1]
}; # end Dstep
# -----------------------------------------------------------------------------------------------------
if {$ok != 0 } {
	puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
} else {
	puts [format $fmt1 "DONE"  $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
}


set clockF [clock milliseconds];
set Etime [expr {($clockF-$clockI)/1000.}]; 	#seconds
puts "*********************************************"
puts "End transient analysis, Elapsed time=$Etime s"
puts "*********************************************"


# exit
