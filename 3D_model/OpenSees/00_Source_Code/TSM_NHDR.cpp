/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.8 $
// $Date: 2007-02-02 01:30:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/NewElement.cpp,v $
                                                                        
// Written: fmk 
// Created: 08/01
//
// Description: This file contains the implementation for the NewElement class.
//
// What: "@(#) NewElement.cpp, revA"


#include "TSM_NHDR.h"

#include <Information.h>
#include <ElementResponse.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <OPS_Globals.h>
#include <OPS_Stream.h>

#include <elementAPI.h>
#include <UniaxialMaterial.h>

#define PI 3.14159l

static int numMyBearing = 0;

Matrix TSM_NHDR::theMatrix(12, 12);
// Vector TSM_NHDR::theVector(6);

void* OPS_TSM_NHDR()
{
	if (numMyBearing == 0) {
		opserr << "TSM_NHDR element - Written by Jose Gallardo, PUC, 2025\n";
		//opserr << "Based on the paper of Gallardo et al. 2025\n";
		numMyBearing++;
	}

	int ndf = OPS_GetNDF();
	if (ndf != 6) {
		opserr << "WARNING invalid ndf: " << ndf;
		opserr << ", for 3D problem need 6 dfs\n";
		return 0;
	}


	if (OPS_GetNumRemainingInputArgs() < 39) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "Want: TSM_NHDR eleTag iNode jNode kvi Kb Kr Kt Do Di Tr -cavitation fc a kappa phi_m fas -shear a1 a2 a3 fs1 ps1 fs2 ps2 fs3 ps3 fm pm ky fy beta eta phi_max P_phi fcp -orient y1 y2 y3  <-mass m> \n";
		return 0;
	}

	// tags
	int idata[3];
	int num = 3;
	if (OPS_GetIntInput(&num, idata) < 0) {
		opserr << "WARNING: invalid integer inputs\n";
		return 0;
	}

	// data1
	double data1[7];
	num = 7;
	if (OPS_GetDoubleInput(&num, data1) < 0) {
		opserr << "WARNING: invalid double inputs\n";
		return 0;
	}

	// data2
	double data2[5];
	const char* type = OPS_GetString();
	if (strcmp(type, "-cavitation") != 0) {
		opserr << "WARNING: want -cavitation\n";
		return 0;
	}
	num = 5;
	if (OPS_GetDoubleInput(&num, data2) < 0) {
		opserr << "WARNING: invalid cavitation parameters\n";
		return 0;
	}
	
	// data3
	double data3[18];
	type = OPS_GetString();
	if (strcmp(type, "-shear") != 0) {
		opserr << "WARNING: want -shear\n";
		return 0;
	}
	num = 18;
	if (OPS_GetDoubleInput(&num, data3) < 0) {
		opserr << "WARNING: invalid shear parameters\n";
		return 0;
	}

	Vector y(3);
	type = OPS_GetString();
	if (strcmp(type, "-orient") == 0) {
		if (OPS_GetNumRemainingInputArgs() < 3) {
			opserr << "WARNING: insufficient arguments after -orient\n";
			return 0;
		}
		num = 3;
		if (OPS_GetDoubleInput(&num, &y(0)) < 0) {
			opserr << "WARNING: invalid orient value\n";
			return 0;
		}
	}

	// options
	double mass = 0.0;

	if (OPS_GetNumRemainingInputArgs() < 1) {
		return new TSM_NHDR(idata[0], idata[1], idata[2],
			data1[0], data1[1], data1[2], data1[3], data1[4], data1[5], data1[6],
			data2[0], data2[1], data2[2], data2[3], data2[4],
			data3[0], data3[1], data3[2], data3[3], data3[4], data3[5], data3[6], data3[7], data3[8],
			data3[9], data3[10], data3[11], data3[12], data3[13], data3[14], data3[15], data3[16], data3[17],
			y, mass);
	}
	
	while (OPS_GetNumRemainingInputArgs() > 0) {
		
		if (strcmp(type, "-mass") == 0) {
			if (OPS_GetNumRemainingInputArgs() < 1) {
				opserr << "WARNING: insufficient args\n";
				return 0;
			}
			num = 1;
			if (OPS_GetDoubleInput(&num, &mass) < 0) {
				opserr << "WARNING: invalid mass\n";
				return 0;
			}
		}
		
	}
	
	return new TSM_NHDR(idata[0], idata[1], idata[2],
		data1[0], data1[1], data1[2], data1[3], data1[4], data1[5], data1[6],
		data2[0], data2[1], data2[2], data2[3], data2[4],
		data3[0], data3[1], data3[2], data3[3], data3[4], data3[5], data3[6], data3[7], data3[8],
		data3[9], data3[10], data3[11], data3[12], data3[13], data3[14], data3[15], data3[16], data3[17],
		y, mass);
}

// constructors:
TSM_NHDR::TSM_NHDR(int tag, int Nd1, int Nd2,
    double _Kvi, double _Kb, double _Kr, double _Kt, double _Do, double _Di,
	double _Tr, double _fc, double _a, double _k, double _pm, double _fas,
	double _a1, double _a2, double _a3, double _fs1, double _ps1,
	double _fs2, double _ps2, double _fs3, double _ps3, double _fm, double _p_m, double _ky,
	double _fy, double _beta, double _eta, double _phimax, double _Pphi, double _fcp,
	const Vector _y, double _mass)
 : Element(tag,ELE_TAG_TSM_NHDR),
  connectedExternalNodes(2), Kvi(_Kvi), fas(_fas), y(_y), Do(_Do), Di(_Di),
	Tr(_Tr), fc(_fc), a(_a), k(_k), pm(_pm),
    Kb(_Kb), Kr(_Kr), Kt(_Kt), 
	a1(_a1), a2(_a2), a3(_a3),
	fs1(_fs1), ps1(_ps1), fs2(_fs2), ps2(_ps2),
	fs3(_fs3), ps3(_ps3), fm(_fm), p_m(_p_m),
	ky(_ky), fy(_fy), beta(_beta), eta(_eta),
	phi_max(_phimax), P_phi(_Pphi), fcp(_fcp),
	L(0.0), ub(6), qb(6), kt(6, 6), mass(_mass),
    ul(12), Tgl(12, 12), Tlb(6, 12), ub_c(6), kinit(12,12), theLoad(12)
{
	
	// ensure the connectedExternalNode ID is of correct size & set values
	if (connectedExternalNodes.Size() != 2) {
		opserr << "TSM_NHDR::TSM_NHDR() - element: "
			<< this->getTag() << " - failed to create an ID of size 2.\n";
		exit(-1);
	}

	connectedExternalNodes(0) = Nd1;
	connectedExternalNodes(1) = Nd2;

	// set node pointers to NULL
	for (int i = 0; i < 2; i++)
		theNodes[i] = 0;

	// initialize initial stiffness matrix
	kinit.Zero();
	kinit.resize(6, 6);
	kinit(0, 0) = Kvi;
	kinit(1, 1) = a1;
	kinit(2, 2) = a1;
	kinit(3, 3) = Kt;
	kinit(4, 4) = Kb;
	kinit(5, 5) = Kb;

	ub.resize(6);
	qb.resize(6);
	ub_c.resize(6);
	
	// initialize other variables
	this->revertToStart();
	
}

TSM_NHDR::TSM_NHDR()
 :Element(0,ELE_TAG_TSM_NHDR),
  connectedExternalNodes(2), Kvi(0.0), y(0), Do(0.0), Di(0.0),
	Tr(0.0), fc(0.0), a(0.0), k(0.0), pm(0.0), fas(0.0),
	Kb(0.0), Kr(0.0), Kt(0.0), 
	a1(0.0), a2(0.0), a3(0.0),
	fs1(0.0), ps1(0.0), fs2(0.0), ps2(0.0),
	fs3(0.0), ps3(0.0), fm(0.0), p_m(0.0),
	ky(0.0), fy(0.0), beta(0.0), eta(0.0),
	phi_max(0.0), P_phi(0.0), fcp(0.0),
	alpha(0.0), phis(0.0), uip(0.0), zx(0.0), zy(0.0),
	t_phis(0.0), t_uip(0.0), t_zx(0.0), t_zy(0.0),
	L(0.0), ub(6), qb(6), kt(6, 6), mass(0.0),
	ul(12), Tgl(12, 12), Tlb(6, 12), ub_c(6), kinit(6, 6), theLoad(12),
	s1(0.0), s2(0.0), t_s1(0.0), t_s2(0.0), theta1(0.0), theta2(0.0), t_theta1(0.0), t_theta2(0.0),
	fpc(0.0), t_fpc(0.0), um(0.0), t_um(0.0)
{
	// ensure the connectedExternalNode ID is of correct size
	if (connectedExternalNodes.Size() != 2) {
		opserr << "TSM_NHDR::TSM_NHDR() - element: "
			<< this->getTag() << " - failed to create an ID of size 2.\n";
		exit(-1);
	}

	// set node pointers to NULL
	for (int i = 0; i < 2; i++)
		theNodes[i] = 0;

}

//  destructor:
TSM_NHDR::~TSM_NHDR()
{

}

int TSM_NHDR::getNumExternalNodes(void) const
{
	//opserr << "getNumExternalNodes \n";
    return 2;
}

const ID & TSM_NHDR::getExternalNodes(void)
{
	//opserr << "getExternalNodes \n";
    return connectedExternalNodes;
}

Node** TSM_NHDR::getNodePtrs(void)
{
	//opserr << "getNodePtrs \n";
  return theNodes;
}

int TSM_NHDR::getNumDOF(void)
{
	//opserr << "getNumDoF \n";
    return 12;
}


void TSM_NHDR::setDomain(Domain *theDomain)
{
	//opserr << "setDomain \n";
	// check Domain is not null - invoked when object removed from a domain
	if (!theDomain) {
		theNodes[0] = 0;
		theNodes[1] = 0;
		return;
	}

	// first set the node pointers
	theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
	theNodes[1] = theDomain->getNode(connectedExternalNodes(1));

	// if can't find both - send a warning message
	if (!theNodes[0] || !theNodes[1]) {
		if (!theNodes[0]) {
			opserr << "WARNING TSM_NHDR::setDomain() - Nd1: "
				<< connectedExternalNodes(0)
				<< " does not exist in the model for";
		}
		else {
			opserr << "WARNING TSM_NHDR::setDomain() - Nd2: "
				<< connectedExternalNodes(1)
				<< " does not exist in the model for";
		}
		opserr << " element: " << this->getTag() << ".\n";

		return;
	}

	// now determine the number of dof and the dimension
	int dofNd1 = theNodes[0]->getNumberDOF();
	int dofNd2 = theNodes[1]->getNumberDOF();

	// if differing dof at the ends - print a warning message
	if (dofNd1 != 6) {
		opserr << "TSM_NHDR::setDomain() - node 1: "
			<< connectedExternalNodes(0)
			<< " has incorrect number of DOF (not 6).\n";
		return;
	}
	if (dofNd2 != 6) {
		opserr << "TSM_NHDR::setDomain() - node 2: "
			<< connectedExternalNodes(1)
			<< " has incorrect number of DOF (not 6).\n";
		return;
	}

	// call the base class method
    this->DomainComponent::setDomain(theDomain);
	// set up the transformation matrix for orientation

	this->setUp();
}   	 

int TSM_NHDR::commitState()
{
	//opserr << "CommitState \n";
	int retVal = 0;

	s1 = t_s1;
	s2 = t_s2;
	theta1 = t_theta1;
	theta2 = t_theta2;
	ub_c = ub;
	um = t_um;
	fpc = t_fpc;



	gamma = t_gamma;
	gammac = t_gammac;
	gammal = t_gammal;

	zx = t_zx;
	zy = t_zy;
	uip = t_uip;
	phis = t_phis;

	return retVal;
}

int TSM_NHDR::revertToLastCommit()
{
	//opserr << "revertToLastCommit \n";
	int retVal = 0;
	
	return retVal;
}

int TSM_NHDR::revertToStart()
{


	double Ks = 25.0;    // Delete 



	//opserr << "revertToStart \n";
	int errCode = 0;

	ub.Zero();
	qb.Zero();
	ub_c.Zero();


	ndir = 360;
	Vector angle(ndir);

	angle(0) = 0.0;
	double inc_angle = 2 * PI / ndir;
	dirx.resize(ndir);
	diry.resize(ndir);
	dirx(0) = 1.0;
	diry(0) = 0.0;

	for (int i = 1; i < ndir; i++) {
		angle(i) = angle(i - 1) + inc_angle;
		dirx(i) = cos(angle(i));
		diry(i) = sin(angle(i));
	}

	gamma.resize(ndir);
	gammac.resize(ndir);
	gammal.resize(ndir);
	t_gamma.resize(ndir);
	t_gammac.resize(ndir);
	t_gammal.resize(ndir);
	gamma.Zero();
	gammac.Zero();
	gammal.Zero();
	t_gamma.Zero();
	t_gammac.Zero();
	t_gammal.Zero();

	zx = 0.0;
	zy = 0.0;
	uip = 0.0;
	phis = 0.0;

	t_zx = 0.0;
	t_zy = 0.0;
	t_uip = 0.0;
	t_phis = 0.0;

	s1 = 0.0;
	s2 = 0.0;
	theta1 = 0.0;
	theta2 = 0.0;
	t_s1 = 0.0;
	t_s2 = 0.0;
	t_theta1 = 0.0;
	t_theta2 = 0.0;
	um = fc / Kvi;
	t_um = um;
	fpc = fc;
	t_fpc = fpc;

	// reset stiffness matrix in basic system
	kt = kinit;

	return errCode;
}

int TSM_NHDR::update(void)
{
	//opserr << "update \n";
	int errCode = 0;
	// get global trial displacements and velocities
	const Vector& dsp1 = theNodes[0]->getTrialDisp();
	const Vector& dsp2 = theNodes[1]->getTrialDisp();
	const Vector& vel1 = theNodes[0]->getTrialVel();
	const Vector& vel2 = theNodes[1]->getTrialVel();

	const Vector& dsp1a = theNodes[0]->getDisp();
	const Vector& dsp2a = theNodes[1]->getDisp();


	static Vector ug(12), ugdot(12), uldot(12), ubdot(6), uga(12), ula(12), uba(6);
	for (int i = 0; i < 6; i++) {
		ug(i) = dsp1(i);  ugdot(i) = vel1(i);
		ug(i + 6) = dsp2(i);  ugdot(i + 6) = vel2(i);
		uga(i) = dsp1(i);
		uga(i + 6) = dsp2a(i);
	}

	// transform response from the global to the local system
	ul.addMatrixVector(0.0, Tgl, ug, 1.0);
	uldot.addMatrixVector(0.0, Tgl, ugdot, 1.0);
	ula.addMatrixVector(0.0, Tgl, uga, 1.0);

	// transform response from the local to the basic system
	ub.addMatrixVector(0.0, Tlb, ul, 1.0);
	ubdot.addMatrixVector(0.0, Tlb, uldot, 1.0);
	uba.addMatrixVector(0.0, Tlb, ula, 1.0);

	// 1) get axial force and stiffness in basic x-direction

	double uh = sqrt(pow(ub(1),2)+ pow(ub(2),2));
	double ac = Di / Do;
	double R = sqrt(1 + pow(ac, 2)) * Do / 4;

	double Kv;
	Kv = Kvi / (1 + (fas / pow(PI, 2)) * pow(uh / R, 2));
	kt.resize(6, 6);
	qb.resize(6);
	qb.Zero();
	
	
	if (ub(0) > 0.0) {
		t_um = max(ub(0), um);
		double uc = fc / Kv;
		double phi = pm * (1 - expo(-ac * ((t_um - uc) / uc)));
		double fcn = fc * (1 - phi);
		double ucn = fcn / Kv;
		double kd = (fpc - fcn) / (t_um - ucn);

		if (ub(0) <= ucn) {
			qb(0) = Kv * ub(0);
			kt(0, 0) = Kv;
		}
		else if (ub(0) < um) {
			qb(0) = fcn + kd * (ub(0) - ucn);
			kt(0, 0) = kd;
		}
		else {
			qb(0) = fc * (1 + (1 / (k * Tr)) * (1 - expo(-k * (ub(0) - uc))));
			kt(0, 0) = (fc/k)* expo(-k * (ub(0) - uc));
			t_fpc = max(fpc, qb(0));
		}
	}
	else {
		qb(0) = ub(0) * Kv;
		kt(0, 0) = Kv;
	}

	// 2) evaluate the shear force and stiffness

	double s1_p = s1;
	double s2_p = s2;

	double theta1_p = theta1;
	double theta2_p = theta2;

	t_s1 = s1;
	t_s2 = s2;
	t_theta1 = theta1;
	t_theta2 = theta2;

	int iter = 0;
	double tol = 1.0e-15;
	double error = 1.0;
	int max_iter = 50;
	double fshear1 = 0;
	double fshear2 = 0;
	double fshear1_p = 0;
	double fshear2_p = 0;

	double F1 = 0.0;
	double F2 = 0.0;

	int flag = 0;
	int count = 0;
	

	double zmx = zx;
	double zmy = zy;
	double gamma_B = 1 - beta;
	double dt = ops_Dt;


	double ks1 = 0.0;
	double ks2 = 0.0;
	double ks3 = 0.0;
	double km = 0.0;
	double gu_n = 0;

	do {
		t_s1 = ub(1) - L * t_theta1;
		t_s2 = ub(2) - L * t_theta2;

		double sxp = (t_s1 - s1) / dt;
		double syp = (t_s2 - s2) / dt;

		// Disipative component

		while (flag == 0) {
			if (zmx == 0) { 
				zmx = 1.0e-15;
			};
			if (zmy == 0) { 
				zmy = 1.0e-15;
			};

			double zp_x = (ky / fy) * (sxp/Tr) - (ky / fy) * zmx * ((beta - phis) * abs((sxp/Tr) * zmx) + (gamma_B - phis) * (sxp/Tr) * zmx + (beta - phis) * abs((syp/Tr) * zmy) + (gamma_B - phis) * (syp/Tr) * zmy) * pow(pow(zmx, 2) + pow(zmy, 2), (eta / 2.0 - 1.0));
			double zp_y = (ky / fy) * (syp/Tr) - (ky / fy) * zmy * ((beta - phis) * abs((sxp/Tr) * zmx) + (gamma_B - phis) * (sxp/Tr) * zmx + (beta - phis) * abs((syp/Tr) * zmy) + (gamma_B - phis) * (syp/Tr) * zmy) * pow(pow(zmx, 2) + pow(zmy, 2), (eta / 2.0 - 1.0));

			double zf_x = zx + dt*zp_x;
			double zf_y = zy + dt*zp_y;

			if (pow(pow(zf_x - t_zx, 2) + pow(zf_y - t_zy, 2), 0.5) < 1.0e-6) {
				flag = 1;
			}
			else {
				t_zx = zf_x;
				t_zy = zf_y;
			}
			zmx = zx + 0.5 * dt*zp_x;
			zmy = zy + 0.5 * dt*zp_y;

			if (count > 100) {
				zmx = t_zx;
				zmy = t_zy;
			}
			if (count > 110) {
				flag = 1;
			}
			count++;
		}
		
		double Rx = fy * t_zx;
		double Ry = fy * t_zy;

		// Hyperelastic component

		double gu_n = pow(t_s1 * t_s1 + t_s2 * t_s2, 0.5) / Tr;
		double gu_a = pow(s1 * s1 + s2 * s2, 0.5) / Tr;

		double um = pow(t_s1 * t_s1 + t_s2 * t_s2, 0.5);

		double unx, uny;

		if (um == 0) {
			unx = 0;
			uny = 0;
		}
		else {
			unx = t_s1 / um;
			uny = t_s2 / um;
		}

		Vector ip(ndir);
		ip.Zero();
		ip(0) = t_s1 * dirx(0) + t_s2 * diry(0);
		int pos = 0;
		double value = ip(0);
		for (int i = 1; i < ndir; i++) {
			ip(i) = t_s1 * dirx(i) + t_s2 * diry(i);

			if (ip(i) > value) {
				pos = i;
				value = ip(i);
			}

		}

		double dg = gu_n - gu_a;

		t_gammal(pos) = max(gammal(pos), gu_n);
		double gamma_trial = gamma(pos) + (9 * sgn(dg) - 7) * (dg / 4.0);

		t_gamma(pos) = min(t_gammal(pos), gamma_trial);
		t_gammac(pos) = gammac(pos) - (1 - sgn(dg)) * (dg / 2.0);

		double ang = 0;

		double km = exp(-fm * pow(gammac(pos), p_m));

		for (int i = 0; i < ndir; i++) {
			ang = dirx(pos) * dirx(i) + diry(pos) * diry(i);
			t_gammal(i) = max((fcp + (1.0 - fcp) * ang) * t_gammal(pos), gammal(i));
			t_gamma(i) = max((fcp + (1.0 - fcp) * ang) * t_gamma(pos), gamma(i));
			t_gammac(i) = max((fcp + (1.0 - fcp) * ang) * t_gammac(pos), gammac(i));
		}

		ks1 = exp(-fs1 * pow(t_gamma(pos), ps1));
		ks2 = ks1 * exp(-fs2 * pow(t_gamma(pos), ps2));
		ks3 = ks1 * exp(-fs3 * pow(t_gamma(pos), ps3));

		km = exp(-fm * pow(t_gammac(pos), p_m));

		double frs1 = unx * (a1 * km * ks1 * gu_n - a2 * ks2 * pow(gu_n, 3.0) + a3 * ks3 * pow(gu_n, 5.0));
		double frs2 = uny * (a1 * km * ks1 * gu_n - a2 * ks2 * pow(gu_n, 3.0) + a3 * ks3 * pow(gu_n, 5.0));

		double fshear1 = Rx + frs1;
		double fshear2 = Ry + frs2;
		
		F1 = fshear1 + qb(0) * t_theta1;
		F2 = fshear2 + qb(0) * t_theta2;

		t_theta1 = (F1 * L - qb(0) * t_s1) / (Kr + qb(0) * L);
		t_theta2 = (F2 * L - qb(0) * t_s2) / (Kr + qb(0) * L);

		error = pow(abs(fshear1_p - fshear1)+ abs(fshear2_p - fshear2) + abs(s1_p - t_s1)+ abs(s2_p - t_s2) + abs(theta1_p - t_theta1)+ abs(theta2_p - t_theta2), 2);
		s1_p = t_s1;
		s2_p = t_s2;
		theta1_p = t_theta1;
		theta2_p = t_theta2;
		fshear1_p = fshear1;
		fshear2_p = fshear2;
		iter++;
	} while ((error > tol) && (iter < max_iter));


	double Nu = pow(t_s1 * t_s1 + t_s2 * t_s2, 0.5);
	double Nua = pow(s1 * s1 + s2 * s2, 0.5);

	if (Nu < Nua) {
		t_uip = max(Nua - fy / ky, uip);
		t_phis = phi_max * (1.0 - exp(-P_phi * abs(ky * t_uip / fy)));
	}


	qb(1) = F1;
	qb(2) = F2;

	double dg_du = 1 / Tr;
	double dRx_dg = ky - ky * t_zx * pow(pow(zmx, 2) + pow(zmy, 2), (eta / 2.0 - 1.0)) * ((beta - phis) * abs(t_zx) + (gamma_B - phis) * t_zx);
	double dRy_dg = ky - ky * t_zy * pow(pow(zmx, 2) + pow(zmy, 2), (eta / 2.0 - 1.0)) * ((beta - phis) * abs(t_zy) + (gamma_B - phis) * t_zy);
	double dfrsx_dg = a1 * km * ks1 - 3.0 * a2 * ks2 * pow(gu_n, 2.0) + 5.0 * a3 * ks3 * pow(gu_n, 4.0);
	double dfrsy_dg = a1 * km * ks1 - 3.0 * a2 * ks2 * pow(gu_n, 2.0) + 5.0 * a3 * ks3 * pow(gu_n, 4.0);

	kt(1, 1) = (dRx_dg + dfrsx_dg) * dg_du - qb(0) / (Kr + qb(0) * L);
	kt(2, 2) = (dRy_dg + dfrsy_dg) * dg_du - qb(0) / (Kr + qb(0) * L);


	// 3) Torsion: force and stiffness
	
	qb(3) = Kt * ub(3);
	kt(3, 3) = Kt;

	// 4) bending: force and stiffness

	qb(4) = ub(4) * Kb;
	qb(5) = ub(5) * Kb;
	kt(4, 4) = Kb;
	kt(5, 5) = Kb;
	
  return 0;
}

const Matrix& TSM_NHDR::getTangentStiff(void)
{
	//opserr << "getTangentStiff \n";
	// zero the matrix
	theMatrix.Zero();

	// transform from basic to local system
	static Matrix kl(12, 12);
	kl.addMatrixTripleProduct(0.0, Tlb, kt, 1.0);

	// // add geometric stiffness to local stiffness
	// double kGeo1 = 0.5 * qb(0);
	// kl(2, 1) -= kGeo1;
	// kl(2, 4) += kGeo1;
	// kl(5, 1) -= kGeo1;
	// kl(5, 4) += kGeo1;
	// double kGeo2 = kGeo1 * 0.5 * L;
	// kl(2, 2) += kGeo2;
	// kl(5, 2) -= kGeo2;
	// double kGeo3 = kGeo1 * (1.0 - 0.5) * L;
	// kl(2, 5) -= kGeo3;
	// kl(5, 5) += kGeo3;

	// transform from local to global system
	static Matrix kg(12, 12);
	kg.addMatrixTripleProduct(0.0, Tgl, kl, 1.0);
	
	return kg;
}

const Matrix& TSM_NHDR::getInitialStiff(void)
{
	//opserr << "getInitialStiff \n";
	// zero the matrix
	theMatrix.Zero();

	// transform from basic to local system
	static Matrix klInit(12, 12);
	klInit.addMatrixTripleProduct(0.0, Tlb, kinit, 1.0);

	// transform from local to global system
	theMatrix.addMatrixTripleProduct(0.0, Tgl, klInit, 1.0);

	return theMatrix;
}
    
void TSM_NHDR::zeroLoad(void)
{
	//opserr << "zeroLoad \n";
	theLoad.Zero();
}

int TSM_NHDR::addLoad(const Vector &addP)
{
	// opserr << "addLoad \n";
	opserr << "TSM_NHDR::addLoad() - "
		<< "load type unknown for element: "
		<< this->getTag() << ".\n";

	return -1;
}

int TSM_NHDR::addInertiaLoadToUnbalance(const Vector &accel)
{
	//opserr << "addInertiaLoadToUnbalance \n";
	// check for quick return
	if (mass == 0.0)
		return 0;

	// get R * accel from the nodes
	const Vector& Raccel1 = theNodes[0]->getRV(accel);
	const Vector& Raccel2 = theNodes[1]->getRV(accel);

	if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
		opserr << "TSM::addInertiaLoadToUnbalance() - "
			<< "matrix and vector sizes are incompatible.\n";
		return -1;
	}

	// want to add ( - fact * M R * accel ) to unbalance
	// take advantage of lumped mass matrix
	double m = 0.5 * mass;
	for (int i = 0; i < 3; i++) {
		theLoad(i) -= m * Raccel1(i);
		theLoad(i + 6) -= m * Raccel2(i);
	}

	return 0;
}

const Vector& TSM_NHDR::getResistingForce()
{	
	//opserr << "getResistingForce\n";
	// zero the residual
	theVector.Zero();
	theVector.resize(12);
	// determine resisting forces in local system
	static Vector ql(12);
	ql.addMatrixTransposeVector(0.0, Tlb, qb, 1.0);

	// add P-Delta moments to local forces
	// double kGeo1 = 0.5 * qb(0);
	// double MpDelta1 = kGeo1 * (ul(4) - ul(1));
	// ql(2) += MpDelta1;
	// ql(5) += MpDelta1;
	// double MpDelta2 = kGeo1 * 0 * L * ul(2);
	// ql(2) += MpDelta2;
	// ql(5) -= MpDelta2;
	// double MpDelta3 = kGeo1 * (1.0 - 0) * L * ul(5);
	// ql(2) -= MpDelta3;
	// ql(5) += MpDelta3;

	// determine resisting forces in global system
	theVector.addMatrixTransposeVector(0.0, Tgl, ql, 1.0);

  return theVector;
}

const Vector& TSM_NHDR::getResistingForceIncInertia()
{
	//opserr << "getResistingForceIncInertia \n";
	// this already includes damping forces from materials
	theVector = this->getResistingForce();

	// subtract external load
	theVector.addVector(1.0, theLoad, -1.0);

	// add inertia forces from element mass
	if (mass != 0.0) {
		const Vector& accel1 = theNodes[0]->getTrialAccel();
		const Vector& accel2 = theNodes[1]->getTrialAccel();

		double m = 0.5 * mass;
		for (int i = 0; i < 3; i++) {
			theVector(i) += m * accel1(i);
			theVector(i + 6) += m * accel2(i);
		}
	}

	return theVector;
}

// modificar ------------------------------------------------------------------------------------------------------------------------------------------------------

int TSM_NHDR::sendSelf(int commitTag, Channel &theChannel)
{
	//opserr << " sendSelf \n";
	static Vector data(17);
	data(0) = this->getTag();
	data(1) = Kvi;
	data(2) = Kb;
	data(3) = Kr;
	data(4) = Kt;
	data(5) = Do;
	data(6) = Di;
	data(7) = Tr;
	data(8) = fc;
	data(9) = a;
	data(10) = k;
	data(11) = pm;
	data(12) = fas;
	data(13) = Ks;
	data(14) = mass;
	data(15) = x.Size();
	data(16) = y.Size();
	theChannel.sendVector(0, commitTag, data);

	// send the two end nodes
	theChannel.sendID(0, commitTag, connectedExternalNodes);

	// send remaining data
	if (x.Size() == 3)
		theChannel.sendVector(0, commitTag, x);
	if (y.Size() == 3)
		theChannel.sendVector(0, commitTag, y);

  return 0;
}

int TSM_NHDR::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{

	// receive element parameters
	static Vector data(17);
	theChannel.recvVector(0, commitTag, data);
	this->setTag((int)data(0));
	Kvi = data(1);
	Kb = data(2);
	Kr = data(3);
	Kt = data(4);
	Do = data(5);
	Di = data(6);
	Tr = data(7);
	fc = data(8);
	a = data(9);
	k = data(10);
	pm = data(11);
	fas = data(12);
	Ks = data(13);
	mass = data(14);

	// receive the two end nodes
	theChannel.recvID(0, commitTag, connectedExternalNodes);


	// receive remaining data
	if ((int)data(15) == 3) {
		x.resize(3);
		theChannel.recvVector(0, commitTag, x);
	}
	if ((int)data(16) == 3) {
		y.resize(3);
		theChannel.recvVector(0, commitTag, y);
	}
	

	// initialize initial stiffness matrix
	kinit.Zero();
	kinit.resize(6, 6);
	kinit(0, 0) = Kvi;
	kinit(1, 1) = Ks;
	kinit(2, 2) = Ks;
	kinit(3, 3) = Kt;
	kinit(4, 4) = Kb;
	kinit(5, 5) = Kb;

	// initialize other variables
	this->revertToStart();

	return 0;
}


int TSM_NHDR::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	//opserr << "displaySelf \n";
	// first determine the end points of the element based on
	// the display factor (a measure of the distorted image)
	const Vector& end1Crd = theNodes[0]->getCrds();
	const Vector& end2Crd = theNodes[1]->getCrds();

	static Vector v1(3);
	static Vector v2(3);

	if (displayMode >= 0) {
		const Vector& end1Disp = theNodes[0]->getDisp();
		const Vector& end2Disp = theNodes[1]->getDisp();

		for (int i = 0; i < 3; i++) {
			v1(i) = end1Crd(i) + end1Disp(i) * fact;
			v2(i) = end2Crd(i) + end2Disp(i) * fact;
		}
	}
	else {
		int mode = displayMode * -1;
		const Matrix& eigen1 = theNodes[0]->getEigenvectors();
		const Matrix& eigen2 = theNodes[1]->getEigenvectors();

		if (eigen1.noCols() >= mode) {
			for (int i = 0; i < 3; i++) {
				v1(i) = end1Crd(i) + eigen1(i, mode - 1) * fact;
				v2(i) = end2Crd(i) + eigen2(i, mode - 1) * fact;
			}
		}
		else {
			for (int i = 0; i < 3; i++) {
				v1(i) = end1Crd(i);
				v2(i) = end2Crd(i);
			}
		}
	}

	return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag(), 0);
}

void TSM_NHDR::Print(OPS_Stream &s, int flag)
{
	//opserr << "Print \n";
	if (flag == OPS_PRINT_CURRENTSTATE) {
		// print everything
		s << "Element: " << this->getTag() << endln;
		s << "  type: TSM\n";
		s << "  iNode: " << connectedExternalNodes(0);
		s << "  jNode: " << connectedExternalNodes(1) << endln;
		//s << "  Shear Material: " << theMaterials[0]->getTag() << endln;
		s << "  mass: " << mass << endln;
		// determine resisting forces in global system
		s << "  resisting force: " << this->getResistingForce() << endln;
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": " << this->getTag() << ", ";
		s << "\"type\": \"TSM\", ";
		s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
		//s << "\"Shear material\": [\"";
		//s << theMaterials[0]->getTag() << "\"], ";
		s << "\"mass\": " << mass << "}";
	}
}

Response* TSM_NHDR::setResponse(const char **argv, int argc, OPS_Stream &output)
{
	//opserr << "setRespone\n";
	Response* theResponse = 0;

	output.tag("ElementOutput");
	output.attr("eleType", "TSM");
	output.attr("eleTag", this->getTag());
	output.attr("node1", connectedExternalNodes[0]);
	output.attr("node2", connectedExternalNodes[1]);

	// global forces
	if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
		strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0)
	{
		output.tag("ResponseType", "Px_1");
		output.tag("ResponseType", "Py_1");
		output.tag("ResponseType", "Pz_1");
		output.tag("ResponseType", "Mx_1");
		output.tag("ResponseType", "My_1");
		output.tag("ResponseType", "Mz_1");
		output.tag("ResponseType", "Px_2");
		output.tag("ResponseType", "Py_2");
		output.tag("ResponseType", "Pz_2");
		output.tag("ResponseType", "Mx_2");
		output.tag("ResponseType", "My_2");
		output.tag("ResponseType", "Mz_2");

		theResponse = new ElementResponse(this, 1, theVector);
	}
	// local forces
	else if (strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "localForces") == 0)
	{
		output.tag("ResponseType", "N_1");
		output.tag("ResponseType", "V1_1");
		output.tag("ResponseType", "V2_2");
		output.tag("ResponseType", "T_1");
		output.tag("ResponseType", "M1_1");
		output.tag("ResponseType", "M1_1");
		output.tag("ResponseType", "N_2");
		output.tag("ResponseType", "V1_2");
		output.tag("ResponseType", "V2_2");
		output.tag("ResponseType", "T_2");
		output.tag("ResponseType", "M1_2");
		output.tag("ResponseType", "M2_2");

		theResponse = new ElementResponse(this, 2, theVector);
	}
	// basic forces
	else if (strcmp(argv[0], "basicForce") == 0 || strcmp(argv[0], "basicForces") == 0)
	{
		output.tag("ResponseType", "qb1");
		output.tag("ResponseType", "qb2");
		output.tag("ResponseType", "qb3");
		output.tag("ResponseType", "qb4");
		output.tag("ResponseType", "qb5");
		output.tag("ResponseType", "qb6");

		theResponse = new ElementResponse(this, 3, Vector(6));
	}
	// local displacements
	else if (strcmp(argv[0], "localDisplacement") == 0 ||
		strcmp(argv[0], "localDisplacements") == 0)
	{
		output.tag("ResponseType", "ux_1");
		output.tag("ResponseType", "uy_1");
		output.tag("ResponseType", "uz_1");
		output.tag("ResponseType", "rx_1");
		output.tag("ResponseType", "ry_1");
		output.tag("ResponseType", "rz_1");
		output.tag("ResponseType", "ux_2");
		output.tag("ResponseType", "uy_2");
		output.tag("ResponseType", "uz_2");
		output.tag("ResponseType", "rx_2");
		output.tag("ResponseType", "ry_2");
		output.tag("ResponseType", "rz_2");

		theResponse = new ElementResponse(this, 4, theVector);
	}
	// basic displacements
	else if (strcmp(argv[0], "deformation") == 0 || strcmp(argv[0], "deformations") == 0 ||
		strcmp(argv[0], "basicDeformation") == 0 || strcmp(argv[0], "basicDeformations") == 0 ||
		strcmp(argv[0], "basicDisplacement") == 0 || strcmp(argv[0], "basicDisplacements") == 0)
	{
		output.tag("ResponseType", "ub1");
		output.tag("ResponseType", "ub2");
		output.tag("ResponseType", "ub3");
		output.tag("ResponseType", "ub4");
		output.tag("ResponseType", "ub5");
		output.tag("ResponseType", "ub6");

		theResponse = new ElementResponse(this, 5, Vector(6));
	}
	
	output.endTag(); // ElementOutput

	return theResponse;
}

int TSM_NHDR::getResponse(int responseID, Information &eleInfo)
{
	//opserr << "getRespone\n";
	// double kGeo1, MpDelta1, MpDelta2, MpDelta3;

	switch (responseID) {
	case 1:  // global forces
		return eleInfo.setVector(this->getResistingForce());

	case 2:  // local forces
		theVector.Zero();
		// determine resisting forces in local system
		theVector.addMatrixTransposeVector(0.0, Tlb, qb, 1.0);
		// add P-Delta moments
		// kGeo1 = 0.5 * qb(0);
		// MpDelta1 = kGeo1 * (ul(4) - ul(1));
		// theVector(2) += MpDelta1;
		// theVector(5) += MpDelta1;
		// MpDelta2 = kGeo1 * 0.5 * L * ul(2);
		// theVector(2) += MpDelta2;
		// theVector(5) -= MpDelta2;
		// MpDelta3 = kGeo1 * (1.0 - 0.5) * L * ul(5);
		// theVector(2) -= MpDelta3;
		// theVector(5) += MpDelta3;

		return eleInfo.setVector(theVector);

	case 3:  // basic forces
		return eleInfo.setVector(qb);

	case 4:  // local displacements
		return eleInfo.setVector(ul);

	case 5:  // basic displacements
		return eleInfo.setVector(ub);

	default:
		return -1;
	}
}

int TSM_NHDR::setParameter(const char **argv, int argc, Parameter &param)
{
	//opserr << "setParameter\n";
  return 1;
}   

int TSM_NHDR::updateParameter(int parameterID, Information &info)
{
	// opserr << "updateParameter\n";
  return -1;
}

void TSM_NHDR::setUp()
{
	//opserr << "setUp\n";
	const Vector& end1Crd = theNodes[0]->getCrds();
	const Vector& end2Crd = theNodes[1]->getCrds();
	Vector xp = end2Crd - end1Crd;

	L = xp.Norm();

	if (L > DBL_EPSILON) {
		if (x.Size() == 0) {
			x.resize(3);
			x = xp;
		}
		else {
			opserr << "WARNING TSM::setUp() - "
				<< "element: " << this->getTag()
				<< " - ignoring nodes and using specified "
				<< "local x vector to determine orientation.\n";
		}
	}
	// check that vectors for orientation are of correct size
	if (x.Size() != 3 || y.Size() != 3) {
		opserr << "TSM::setUp() - "
			<< "element: " << this->getTag() << endln
			<< " - incorrect dimension of orientation vectors.\n";
		exit(-1);
	}

	// establish orientation of element for the transformation matrix
	// z = x cross y
	static Vector z(3);
	z(0) = x(1) * y(2) - x(2) * y(1);
	z(1) = x(2) * y(0) - x(0) * y(2);
	z(2) = x(0) * y(1) - x(1) * y(0);

	// y = z cross x
	y(0) = z(1) * x(2) - z(2) * x(1);
	y(1) = z(2) * x(0) - z(0) * x(2);
	y(2) = z(0) * x(1) - z(1) * x(0);

	// compute length(norm) of vectors
	double xn = x.Norm();
	double yn = y.Norm();
	double zn = z.Norm();

	// check valid x and y vectors, i.e. not parallel and of zero length
	if (xn == 0 || yn == 0 || zn == 0) {
		opserr << "TSM::setUp() - "
			<< "element: " << this->getTag() << endln
			<< " - invalid orientation vectors.\n";
		exit(-1);
	}

	// create transformation matrix from global to local system
	Tgl.Zero();
	Tgl(0, 0) = Tgl(3, 3) = Tgl(6, 6) = Tgl(9, 9) = x(0) / xn;
	Tgl(0, 1) = Tgl(3, 4) = Tgl(6, 7) = Tgl(9, 10) = x(1) / xn;
	Tgl(0, 2) = Tgl(3, 5) = Tgl(6, 8) = Tgl(9, 11) = x(2) / xn;
	Tgl(1, 0) = Tgl(4, 3) = Tgl(7, 6) = Tgl(10, 9) = y(0) / yn;
	Tgl(1, 1) = Tgl(4, 4) = Tgl(7, 7) = Tgl(10, 10) = y(1) / yn;
	Tgl(1, 2) = Tgl(4, 5) = Tgl(7, 8) = Tgl(10, 11) = y(2) / yn;
	Tgl(2, 0) = Tgl(5, 3) = Tgl(8, 6) = Tgl(11, 9) = z(0) / zn;
	Tgl(2, 1) = Tgl(5, 4) = Tgl(8, 7) = Tgl(11, 10) = z(1) / zn;
	Tgl(2, 2) = Tgl(5, 5) = Tgl(8, 8) = Tgl(11, 11) = z(2) / zn;

	// create transformation matrix from local to basic system (linear)
	Tlb.Zero();
	Tlb(0, 0) = Tlb(1, 1) = Tlb(2, 2) = Tlb(3, 3) = Tlb(4, 4) = Tlb(5, 5) = -1.0;
	Tlb(0, 6) = Tlb(1, 7) = Tlb(2, 8) = Tlb(3, 9) = Tlb(4, 10) = Tlb(5, 11) = 1.0;
	Tlb(1, 5) = -0.5 * L;
	Tlb(1, 11) = -(1.0 - 0.5) * L;
	Tlb(2, 4) = -Tlb(1, 5);
	Tlb(2, 10) = -Tlb(1, 11);

}

double TSM_NHDR::abs(double x)
{
	if (x < 0) return -x;
	else return x;
}

double TSM_NHDR::max(double x,double y) {
	if (x < y) return y;
	else return x;
}

double TSM_NHDR::min(double x, double y) {
	if (x < y) return x;
	else return y;
}

double TSM_NHDR::expo(double x) {
	double e = 2.718281828459046;
	return pow(e, x);
}

double TSM_NHDR::sgn(double x) {
	if (x > 0)
		return 1.0;
	else if (x < 0)
		return -1.0;
	else
		return 0.0;
}

const Matrix& TSM_NHDR::getDamp()
{
	// opserr << "getDamp \n";
	// zero the matrix
	theMatrix.Zero();

	// call base class to setup Rayleigh damping
	double factThis = 0.0;


	// now add damping tangent from materials
	static Matrix cb(6, 6);
	cb.Zero();
	cb(0, 0) = 0.0;
	cb(1, 1) = 0.0;
	cb(2, 2) = 0.0;
	cb(3, 3) = 0.0;
	cb(4, 4) = 0.0;
	cb(5, 5) = 0.0;

	// transform from basic to local system
	static Matrix cl(12, 12);
	cl.addMatrixTripleProduct(0.0, Tlb, cb, 1.0);

	// transform from local to global system and add to cg
	theMatrix.addMatrixTripleProduct(factThis, Tgl, cl, 1.0);

	return theMatrix;
}

const Matrix& TSM_NHDR::getMass()
{
	// zero the matrix
	theMatrix.Zero();

	// check for quick return
	if (mass == 0.0) {
		return theMatrix;
	}

	double m = 0.5 * mass;
	for (int i = 0; i < 3; i++) {
		theMatrix(i, i) = m;
		theMatrix(i + 6, i + 6) = m;
	}

	return theMatrix;
}



