# Las unidades son en el SI

# Iniciar
wipe
wipeAnalysis
file mkdir Results;

# -------------------------------------------
#			Generación del modelo
# -------------------------------------------

model basic -ndm 3 -ndf 6 ;#BasicBuilder

# Load of parameters
source Data/Data_TSM.tcl

# Element height
set h 34.3

# Creation of the element
node 1  0.0 0.0 0.0;
node 2  0.0 0.0 $h

# Restraints
fix  1   1 1 1 1 1 1
fix  2   0 0 0 1 1 1

# Axial properties
set fas 3.0
set Kb 461998.9196
set Kr 461998.9196
set Kvi 2438.284
set Kt 461998.9196
set Ks 0.90;
set Do 70.0
set Di 10.0

# Parameters for cavitation
set Tr 17.60
set a 1.0
set k  20.0
set pm2 0.75
set fc 56.55
set bet  0.90;
set eta  0.32;

# Properties of the element
uniaxialMaterial Elastic 5 $Kvi
uniaxialMaterial Elastic 6 $Kb
element  TSM_NHDR 1 1 2 $Kvi  $Kb $Kr $Kt $Do $Di $Tr -cavitation $fc $a $k $pm2 $fas -shear $a1 $a2 $a3 $fs1 $ps1 $fs2 $ps2 $fs3 $ps3 $fm $p_m $ky $fy $bet $eta $pmx $Pfi 0.5 -orient 1.0 0.0 0.0

puts "Model ok"

constraints Transformation
numberer RCM
system BandGeneral
test EnergyIncr 1e-7 25 0;
algorithm Newton;
integrator Newmark 0.5 0.25 ;
analysis Transient;

# Información a guardar
recorder Element -file ./Results/Displ_TSM.txt -time -closeOnWrite -ele 1 localDisplacement
recorder Element -file ./Results/Force_TSM.txt -time -closeOnWrite -ele 1 localForce

set dt 0.001;

timeSeries Path 3 -filePath Data/Displ_hist_D1.tcl -dt $dt -factor 0.7071;
timeSeries Path 4 -filePath Data/Vel_hist_D1.tcl -dt $dt;
timeSeries Path 50 -filePath Data/Force_hist_P1.tcl -dt $dt;

pattern MultipleSupport  20  {
     groundMotion 101 Series -disp 3 -vel 4;
     imposedSupportMotion 2 1 101;
	 imposedSupportMotion 2 2 101;
}

# pattern Plain 21 50 {
# 	load 2 0.0 0.0 -1.0 0.0 0.0 0.0
# }
timeSeries Constant 7 -factor 1
set pa [expr 631.0*0.8]
pattern Plain 1 7 {
	load 2 0.0 0.0 [expr -$pa] 0.0 0.0 0.0
}

set Nsteps 11750;

for {set i 0} {$i<$Nsteps} {incr i} {
	#--------------------
	#Step and time
	#--------------------
	set step [expr {$i+1}];         #step number
    set time [expr {$dt*$step}];    #time (final of each step)
			
	algorithm ModifiedNewton;
	set ok [analyze 1 $dt]

	if {$ok != 0} {
		puts "Trying Modified Newton"
		algorithm ModifiedNewton;
		set ok [analyze 1 $dt]
	}
	if {$ok != 0} {
		puts "Trying Broyden"
		algorithm Broyden 8
		set ok [analyze 1 $dt]
	}
	if {$ok != 0} {
		puts "Trying NewtonWithLineSearch"
		algorithm NewtonLineSearch 0.8 
		set ok [analyze 1 $dt]
	}
	if {$ok != 0} {
		puts "CONVERGENCE FAILED!"
	}
			

	#--------------------
	#Print step in tcl and .txt
	#--------------------
	
	set res [expr $step%100]
	
	if {$res == 0} {
		set t1 [clock milliseconds];		#final time step	
		puts ">>> Converged Step=$step, Time=$time s"
	} 
}


puts "*********************************************"
puts "End transient analysis"
puts "*********************************************"
exit