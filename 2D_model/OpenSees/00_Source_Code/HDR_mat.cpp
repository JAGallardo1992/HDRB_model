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

//
// Description: This file contains the class implementation for 
// NewUniaxialMaterial.

#include <HDR_mat.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>
#include <Matrix.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>
#include <elementAPI.h>

void* OPS_HDR_mat()
{
	int numdata = OPS_GetNumRemainingInputArgs();
	if (numdata < 19) {
		opserr << "WARNING: Insufficient arguments\n";
		return 0;
	}

	int tag;
	numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
		opserr << "WARNING invalid tag\n";
		return 0;
	}

	double data[18] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0 };
	numdata = OPS_GetNumRemainingInputArgs();
	if (numdata > 19) {
		numdata = 19;
	}
	if (OPS_GetDoubleInput(&numdata, data)) {
		opserr << "WARNING invalid double inputs\n";
		return 0;
	}

	UniaxialMaterial* mat = new HDR_mat(tag, data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8], data[9], data[10],
		data[11], data[12], data[13], data[14], data[15], data[16], data[17]);
	if (mat == 0) {
		opserr << "WARNING: failed to create HDR_mat material\n";
		return 0;
	}
	return mat;
}

HDR_mat::HDR_mat(int tag, double v_a1, double v_a2, double v_a3, double v_fs1, double v_ps1,
	double v_fs2, double v_ps2, double v_fs3, double v_ps3, double v_fm, double v_pm, double v_ky,
	double v_fy, double v_beta, double v_eta, double v_phimax, double v_Pphi, double v_h)
	: UniaxialMaterial(tag, MAT_TAG_HDR_mat),
	a1(v_a1), a2(v_a2), a3(v_a3),
	fs1(v_fs1), ps1(v_ps1), fs2(v_fs2), ps2(v_ps2),
	fs3(v_fs3), ps3(v_ps3), fm(v_fm), pm(v_pm),
	ky(v_ky), fy(v_fy), beta(v_beta), eta(v_eta),
	phi_max(v_phimax), P_phi(v_Pphi), h(v_h)
{
	this->revertToStart();

}



HDR_mat::HDR_mat()
	:UniaxialMaterial(0, MAT_TAG_HDR_mat)
{
}

HDR_mat::~HDR_mat()
{

}

int
HDR_mat::setTrialStrain(double strain, double strainRate)
{

	double d_strain;
	double gcum;
	double glim;
	double gamma;
	double ks1;
	double ks2;
	double ks3;
	double km;
	double fh;
	double zp = z;
	double zf = z;
	double zm = z;
	int flag = 0;
	int count = 0;
	double Ro = fy * z;
	double R;
	double uipn;
	double phin;

	t_gp = gp;
	t_gn = gn;
	t_gcp = gcp;
	t_gcn = gcn;
	t_glp = glp;
	t_gln = gln;
	t_z = z;
	t_uy = uy;
	t_ui = ui;
	t_uip = uip;
	t_phi = phi;

	// set the trial strain

	tStrain = strain/h;

	d_strain = tStrain - Strain/h;

	if (tStrain >= 0) {
		gcum = gcp;
		if (d_strain >= 0) {
			glim = max(glp, tStrain);
			gamma = min(gp + 0.5 * d_strain, glim);
		}
		else {
			glim = glp;
			gamma = min(gp - 4.0 * d_strain, glim);
			gcum = gcum - d_strain;
		}
		t_gp = gamma;
		t_gcp = gcum;
		t_glp = glim;
	}
	else {
		gcum = gcn;
		if (d_strain < 0) {
			glim = max(gln, -tStrain);
			gamma = min(glim, gn - 0.5 * d_strain);
		}
		else {
			glim = gln;
			gamma = min(glim, gn + 4 * d_strain);
			gcum = gcum + d_strain;
		}
		t_gn = gamma;
		t_gcn = gcum;
		t_gln = glim;
	};

	ks1 = exp(-fs1 * pow(gamma, ps1));
	ks2 = ks1 * exp(-fs2 * pow(gamma, ps2));
	ks3 = ks1 * exp(-fs3 * pow(gamma, ps3));
	km = exp(-fm * pow(gcum, pm));

	fh = ks1 * km * a1 * tStrain + ks2 * a2 * pow(tStrain, 3.0) + ks3 * a3 * pow(tStrain, 5.0);

	t_z = z;
	
	
	
	
	while (flag == 0) {
		zp = (ky / fy) * d_strain * (1 - pow(abs(zm), eta) * (beta * sgn(d_strain * zm) + alpha - phi * sgn(d_strain) * (sgn(zm) + sgn(d_strain))));
		zf = z + zp;
		if (abs(zf - t_z) < 1.0e-6) {
			flag = 1;
		}
		else {
			t_z = zf;
		}
		zm = z + 0.5 * zp;
		if (count > 100) {
			zm = z;
		}
		if (count > 110) {
			flag = 1;
		}
		count++;

	}

	R = fy * t_z;

	//if (abs(sgn(d_strain) - sgn(d_a)) == 2) {
	//	uipn = abs(tStrain - uy - d_strain);
	//	t_uip = max(uipn, uip);
	//	phin = phi_max * (1 - exp(-P_phi * abs(ky * t_uip / fy)));
	//	t_phi = max(phi, phin);
	//	t_uy = tStrain - R / ky + fy / (ky * (1 - 4 * t_phi)) * sgn(d_strain);
	//}

	//if (abs(sgn(t_z) - sgn(z))) {
	//	t_uy = tStrain + R / ky - fy / (ky * (1 - 4 * phi));
	//	t_ui = tStrain - R / (R - Ro) * d_strain;
	//}
	
	if (abs(sgn(d_strain) - sgn(d_a)) == 2) {
		t_uip = abs(tStrain - uy - d_strain);
	}

	if (abs(sgn(t_z) - sgn(z))) {
		t_uy = uy;
		t_phi = phi_max * (1 - exp(-P_phi * abs(t_uip / uy)));
	}
	
	

	tStress = R + fh;

	double dFh_dg = ks1 * km * a1 - 3.0 * ks2 * a2 * pow(gamma, 2.0) + 5.0 * ks3 * a3 * pow(gamma, 4.0);
	double dg_du = 1 / h;
	double dFd_dg = ky - ky*pow(abs(t_z), eta) * (beta * sgn(d_strain * t_z) + alpha - phi * sgn(d_strain) * (sgn(t_z) + sgn(d_strain)));

	tTangent = dFh_dg * dg_du + dFd_dg * dg_du;

	tStrain = tStrain * h;
	return 0;
}

double
HDR_mat::getStress(void)
{
	return tStress;
}

double
HDR_mat::getTangent(void)
{
	return tTangent;
}

double
HDR_mat::getInitialTangent(void)
{
	// return the initial tangent
	return a1 + ky;
}

double
HDR_mat::getStrain(void)
{
	return tStrain;
}

int
HDR_mat::commitState(void)
{
	d_a = tStrain - Strain;
	Strain = tStrain;
	Stress = tStress;
	Tangent = tTangent;

	// Hystory variabes
	gp = t_gp;
	gn = t_gn;
	gcp = t_gcp;
	gcn = t_gcn;
	glp = t_glp;
	gln = t_gln;
	z = t_z;
	uy = t_uy;
	ui = t_ui;
	uip = t_uip;
	phi = t_phi;

	return 0;
}

int
HDR_mat::revertToLastCommit(void)
{
	tStrain = Strain;
	tStress = Stress;
	tTangent = Tangent;

	t_gp = gp;
	t_gn = gn;
	t_gcp = gcp;
	t_gcn = gcn;
	t_glp = glp;
	t_gln = gln;
	t_z = z;
	t_uy = uy;
	t_ui = ui;
	t_uip = uip;
	t_phi = phi;
	d_a = d_a;

	return 0;
}

int
HDR_mat::revertToStart(void)
{
	tStrain = 0.0;
	tStress = 0.0;
	tTangent = 0.0;
	Strain = 0.0;
	Stress = 0.0;
	Tangent = 0.0;

	alpha = 1.0 - beta;
	// history variables
	gp = 0.0;
	gn = 0.0;
	gcp = 0.0;
	gcn = 0.0;
	glp = 0.0;
	gln = 0.0;
	z = 0.0;
	uy = fy / ky;
	ui = 0.0;
	uip = 0.0;
	phi = 0.0;
	d_a = 0.0;

	t_gp = 0.0;
	t_gn = 0.0;
	t_gcp = 0.0;
	t_gcn = 0.0;
	t_glp = 0.0;
	t_gln = 0.0;
	t_z = 0.0;
	t_uy = fy / ky;
	t_ui = 0.0;
	t_uip = 0.0;
	t_phi = 0.0;

	return 0;
}

UniaxialMaterial*
HDR_mat::getCopy(void)
{
	HDR_mat* theCopy = new HDR_mat(this->getTag(),
		a1, a2, a3, fs1, ps1, fs2, ps2, fs3, ps3, fm, pm,
		ky, fy, beta, eta, phi_max, P_phi, h);

	theCopy->tStrain = tStrain;
	theCopy->Strain = Strain;
	theCopy->tStress = tStress;
	theCopy->tTangent = tTangent;
	theCopy->Stress = Stress;
	theCopy->Tangent = Tangent;

	theCopy->alpha = alpha;

	// history variables
	theCopy->gp = gp;
	theCopy->gn = gn;
	theCopy->gcp = gcp;
	theCopy->gcn = gcn;
	theCopy->glp = glp;
	theCopy->gln = gln;
	theCopy->z = z;
	theCopy->uy = uy;
	theCopy->ui = ui;
	theCopy->uip = uip;
	theCopy->phi = phi;
	theCopy->d_a = d_a;

	theCopy->t_gp = t_gp;
	theCopy->t_gn = t_gn;
	theCopy->t_gcp = t_gcp;
	theCopy->t_gcn = t_gcn;
	theCopy->t_glp = t_glp;
	theCopy->t_gln = t_gln;
	theCopy->t_z = t_z;
	theCopy->t_uy = t_uy;
	theCopy->t_ui = t_ui;
	theCopy->t_uip = t_uip;
	theCopy->t_phi = t_phi;

	return theCopy;
}

int HDR_mat::sendSelf(int cTag, Channel& theChannel)
{
	static Vector data(32);
	data(0) = this->getTag();
	data(1) = a1;
	data(2) = a2;
	data(3) = a3;
	data(4) = fs1;
	data(5) = ps1;
	data(6) = fs2;
	data(7) = ps2;
	data(8) = fs3;
	data(9) = ps3;
	data(10) = fm;
	data(11) = pm;
	data(12) = uy;
	data(13) = fy;
	data(14) = ky;
	data(15) = beta;
	data(16) = alpha;
	data(17) = eta;
	data(18) = phi_max;
	data(19) = P_phi;
	data(20) = gp;
	data(21) = gn;
	data(22) = gcp;
	data(23) = gcn;
	data(24) = glp;
	data(25) = gln;
	data(26) = z;
	data(27) = ui;
	data(28) = uip;
	data(29) = phi;
	data(30) = d_a;
	data(31) = h;

	if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
		opserr << "HDR_mat::sendSelf() - failed to send Vector\n";
		return -1;
	}

	return 0;
}


int HDR_mat::recvSelf(int cTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	static Vector data(32);
	if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
		opserr << "HDR_mat::recvSelf() - failed to recvSelf\n";
		return -1;
	}

	this->setTag((int)data(0));
	a1 = data(1);
	a2 = data(2);
	a3 = data(3);
	fs1 = data(4);
	ps1 = data(5);
	fs2 = data(6);
	ps2 = data(7);
	fs3 = data(8);
	ps3 = data(9);
	fm = data(10);
	pm = data(11);
	uy = data(12);
	fy = data(13);
	ky = data(14);
	beta = data(15);
	alpha = data(16);
	eta = data(17);
	phi_max = data(18);
	P_phi = data(19);
	gp = data(20);
	gn = data(21);
	gcp = data(22);
	gcn = data(23);
	glp = data(24);
	gln = data(25);
	z = data(26);
	ui = data(27);
	uip = data(28);
	phi = data(29);
	d_a = data(31);
	h = data(32);
	return 0;
}

void HDR_mat::Print(OPS_Stream& s, int flag)
{
	s << "HDR material with anisotropic damage, tag : " << this->getTag() << endln;
	s << "Developed by J. Gallardo and J. C. de la Llera" << endln;
	return;
}

double HDR_mat::sgn(double x)
{
	if (x > 0)
		return 1.0;
	else if (x < 0)
		return -1.0;
	else
		return 0.0;
}

double HDR_mat::abs(double x)
{
	if (x < 0) return -x;
	else return x;
}

double HDR_mat::max(double x, double y)
{
	if (x >= y) return x;
	else return y;
}

double HDR_mat::min(double x, double y)
{
	if (x >= y) return y;
	else return x;
}

double HDR_mat::exp(double x)
{
	double e = 2.718281828459046;
	return pow(e, x);
}