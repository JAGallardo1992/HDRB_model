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

// $Revision: 1.6 $
// $Date: 2022-01-01 00:41:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/NewHDR.h,v $

//
// Written: Jose Gallardo (jogallardo@uc.cl)
//


#ifndef HDR_mat_h
#define HDR_mat_h

#include <UniaxialMaterial.h>
#include <Matrix.h>

class HDR_mat : public UniaxialMaterial
{
public:
	HDR_mat(int tag,
		double a1,
		double a2,
		double a3,
		double fs1,
		double ps1,
		double fs2,
		double ps2,
		double fs3,
		double ps3,
		double fm,
		double pm,
		double ky,
		double fy,
		double beta,
		double eta,
		double phi_max,
		double P_phi,
		double h);
	HDR_mat();
	~HDR_mat();

	const char* getClassType(void) const { return "HDR_mat"; };

	int setTrialStrain(double strain, double strainRate = 0.0);
	double getStrain(void);
	double getStress(void);
	double getTangent(void);
	double getInitialTangent(void);

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	UniaxialMaterial* getCopy(void);

	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel,
		FEM_ObjectBroker& theBroker);

	void Print(OPS_Stream& s, int flag = 0);

protected:

private:
	// private methods
	double sgn(double);
	double abs(double);
	double max(double, double);
	double min(double, double);
	double exp(double);
	double a1;
	double a2;
	double a3;
	double fs1;
	double ps1;
	double fs2;
	double ps2;
	double fs3;
	double ps3;
	double fm;
	double pm;
	double ky;
	double fy;
	double beta;
	double eta;
	double phi_max;
	double P_phi;
	double h;

	double alpha;
	double d_a;
	// history variables
	double gp;
	double gn;
	double gcp;
	double gcn;
	double glp;
	double gln;
	double z;
	double uy;
	double ui;
	double uip;
	double phi;

	double t_gp;
	double t_gn;
	double t_gcp;
	double t_gcn;
	double t_glp;
	double t_gln;
	double t_z;
	double t_uy;
	double t_ui;
	double t_uip;
	double t_phi;

	// Trial and comitted
	double tStrain;
	double tStress;
	double tTangent;
	double Strain;
	double Stress;
	double Tangent;
};

#endif
