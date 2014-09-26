/*
 *  /src/AuxBioCro.c by Fernando Ezequiel Miguez  Copyright (C) 2007 - 2012
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 or 3 of the License
 *  (at your option).
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  A copy of the GNU General Public License is available at
 *  http://www.r-project.org/Licenses/
 *
 *  Part of the code here (sunML, EvapoTrans, SoilEvapo, TempTo and
 *  the *prof functions) are based on code in WIMOVAC. WIMOVAC is
 *  copyright of Stephen Long and Stephen Humphries.
 *  Documentation for WIMOVAC can be found at
 *  http://www.life.illinois.edu/plantbio/wimovac/ (checked 02-13-2010)
 *
 */


/* This file will contain functions which are common to several */
/* routines in the BioCro package. These are functions needed */
/* internally. The normal user will not need them */

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "c4photo.h"
#include "AuxBioCro.h"


/* lightME function. Light Macro Environment */
void lightME(double lat, int DOY, int td)
{

	extern double tmp1[];
	double *ip1;
	ip1 = &tmp1[0];
	double omega, delta0, delta, deltaR;
	double tf, SSin, CCos, PPo;
	double CosZenithAngle, CosHour;
	double CosHourDeg;
	double Idir, Idiff, propIdir, propIdiff;
	const double DTR = M_PI/180;
	const double tsn = 12.0;
	const double alpha = 0.85;
	const double SolarConstant = 2650;
	const double atmP = 1e5;

	omega = lat * DTR;
	delta0 = 360.0 * ((DOY + 10)/365.0);
	delta = -23.5 * cos(delta0*DTR);
	deltaR = delta * DTR;

	tf = (15.0*(td - tsn))*DTR;
	SSin = sin(deltaR) * sin(omega);
	CCos = cos(deltaR) * cos(omega);

	CosZenithAngle = SSin + CCos * cos(tf);
	if(CosZenithAngle < pow(10,-10))
		CosZenithAngle = pow(10,-10);

	CosHour = -tan(omega) * tan(deltaR);
	CosHourDeg = (1/DTR)*CosHour;
	if(CosHourDeg < -57)
		CosHour = -0.994;

	PPo = 1e5 / atmP;
	Idir = SolarConstant * (pow(alpha,(PPo/CosZenithAngle)));
	Idiff = 0.3 * SolarConstant *(1 - pow(alpha,(PPo/CosZenithAngle))) * CosZenithAngle ;

	propIdir = Idir / (Idir + Idiff);
	propIdiff = Idiff / (Idir + Idiff);

	*ip1 = propIdir;
	*(ip1+1) = propIdiff;
	*(ip1+2) = CosZenithAngle;
	*(ip1+3) = Idir;
	*(ip1+4) = Idiff;
}

/* This is the original sunML function written by Fernando Miguez */
/* void sunML(double Idir, double Idiff, double LAI, int nlayers, double cosTheta, double kd, double chil, double heightf) */
/* { */
/* 	extern int sp1, sp2, sp3, sp4, sp5, sp6; */
/* 	extern double layIdir[], layIdiff[], layItotal[], layFsun[], layFshade[], layHeight[]; */
/* 	int i; */
/* 	double k0, k1, k; */
/* 	double LAIi, CumLAI; */
/* 	double Isolar, Idiffuse, Itotal; */
/* 	double Ls,Ls2, Ld, Ld2; */
/* 	double Fsun, Fshade; */

/* 	k0 = cosTheta * sqrt(pow(chil ,2) + pow(tan(acos(cosTheta)),2)); */
/* 	k1 = chil + 1.744*pow((chil+1.183),-0.733); */
/* 	k = k0/k1; */
/* 	if(k<0) */
/* 		k = -k; */

/* 	LAIi = LAI / nlayers; */

/* 	for(i=0;i<nlayers;i++) */
/* 	{ */
/* 		CumLAI = LAIi * (i+1); */
/* 		if(i == 0){ */
/* 			Isolar = Idir; */
/* 			Idiffuse = Idiff; */
/* 		} */
/* 		else{ */
/* 			Isolar = Idir * exp(-k * CumLAI); */
/* 			Idiffuse = Idiff * exp(-kd * CumLAI); */
/* 		} */
/* 		Itotal = Isolar + Idiffuse; */

/* 		Ls = (1-exp(-k*CumLAI))/k; */
/* 		Ld = CumLAI - Ls; */
/* 		Ls2 = (cosTheta * (1-exp(-k*Ls/cosTheta)))/k; */
/* 		Ld2 = CumLAI - Ls2; */
/* 		Fsun = Ls2 /  CumLAI; */
/* 		Fshade = Ld2 / CumLAI; */

/* 		/\* collecting the results *\/ */
/* 		layIdir[sp1++] = Isolar; */
/* 		layIdiff[sp2++] = Idiffuse; */
/* 		layItotal[sp3++] = Itotal; */
/* 		layFsun[sp4++] = Fsun; */
/* 		layFshade[sp5++] = Fshade; */
/* 		layHeight[sp6++] = CumLAI/heightf; */
/* 	} */
/* } */

/* sunML function. Sun-Shade Multi-Layer model */
/* This function was written by Deepak Jaiswal */
/* June 1 2012 */
/* void sunML(double Idir, double Idiff, double LAI, int nlayers, */
/* 	   double cosTheta, double kd, double chil, double heightf) */
/* { */
/* 	extern int sp1, sp2, sp3, sp4, sp5, sp6; */
/* 	extern double layIdir[], layIdiff[], layItotal[], layFsun[], layFshade[], layHeight[]; */
/* 	int i; */
/* 	double k0, k1, k; */
/* 	double LAIi, CumLAI; */
/* 	double Isolar, Idiffuse, Itotal,Iscat,alphascatter; */
/* 	double Ls,Ls2, Ld, Ld2; */
/* 	double Fsun, Fshade,Fcanopy; */
/* 	double Lssaved,Ldsaved; */

/* 	/\*extinction coefficient evaluated here is defnied as ratio of */
/* 	  shadow on the horizontal area projected from a zenith angle Theta. */
/* 	  Earlier version of the code was based on ratio of area projected in */
/* 	  the direction of incoming beam radiation. */
/* 	  Therefore, cosTheta term in numerator is removed */
/* 	*\/ */
/* 	/\** wishlist**\/ */
/* 	/\* currently a default value of kd is used, this value can be calculated, see page 254 of Campbell and Norman *\/ */
/* 	/\* Long Wave radiation balanc has not been included in the canopy environment *\/ */

/* 	alphascatter=0.8; /\* This value is for PAR, for long wave radiation a shorter value = 0.2 is recommended see page  255 of campbell and Norman *\/ */
/* 	k0 =  sqrt(pow(chil ,2) + pow(tan(acos(cosTheta)),2)); */
/* 	k1 = chil + 1.744*pow((chil+1.183),-0.733); */
/* 	k = k0/k1; */
/* 	if(k<0) k = -k; */
/* 	if(k > 1) k = 1; */

/* 	LAIi = LAI / nlayers; */

/* 	for(i=0;i<nlayers;i++) */
/* 	{ */
/* 		CumLAI = LAIi * (i+1/2); */
/* 		if(i == 0) */
/* 		      { */
		
/* 			Iscat=Idir * cosTheta * exp(-k *sqrt(alphascatter)* CumLAI)-Idir * cosTheta * exp(-k * CumLAI); /\* scattered radiation Eq 15.19 from Campbell and Norman pp 259 *\/ */
/* 			Idiffuse = Idiff * exp(-kd * CumLAI) + Iscat; /\* Eq 15.19 from Campbell and Norman pp 259 *\/ */
/* 			Isolar = Idir * k * cosTheta; /\* Avergae direct radiation on sunlit leaves *\/ */
/* 			Itotal = Isolar + Idiffuse; /\* Total radiation *\/ */
	                
/* 			Fcanopy=CumLAI; /\* leaf area of canopy currenly under consideration *\/ */
/* 			/\* old version Ls = (1-exp(-k*CumLAI))/k; /\* To evaluate sunlit leaf area in Fcanopy*\/ */
/*                         Ls = (1-exp(-k*LAIi))/k * exp(-k*CumLAI); */
/* 			Lssaved=Ls;              /\* To use for evaluating sunlit leaves in the lower layer *\/ */
/* 			Ld=LAIi-Ls;          /\* shaded leaf area *\/ */
/* 			Ldsaved=Ld;              /\* To use for evaluating shaded leaves in the lower layer*\/ */
/* 			Fsun=Ls/(Ls+Ld);         /\*fracction of sunlit leaves in the current layer *\/ */
/* 			Fshade=Ld/(Ls+Ld);       /\* fraction of shaded leaves in the current layer *\/ */
/*                         Itotal = Fsun*Isolar + Idiffuse; */
/* 			} */
		
/* 		else */
/* 		{ */
/* 			Iscat=Idir * cosTheta * exp(-k *sqrt(alphascatter)* CumLAI)-Idir * cosTheta * exp(-k * CumLAI); /\* scattered radiation Eq 15.19 from Campbell and Norman pp 259 *\/ */
/* 			Idiffuse = Idiff * exp(-kd * CumLAI) + Iscat; /\* Eq 15.19 from Campbell and Norman pp 259 *\/ */
/* 			Isolar = Idir * k * cosTheta; /\* Avergae direct radiation on sunlit leaves *\/ */
/* 			Itotal = Isolar + Idiffuse; */
	                
/* 			Fcanopy=CumLAI;  /\* leaf area of canopy currenly under consideration *\/ */
/* 			Ls = (1-exp(-k*LAIi))/k * exp(-k*CumLAI); /\* To evaluate  sunlit leaf area  in Fcanopy*\/ */
/* 			Ld=LAIi-Ls; /\* evaluate shaded leaf area in Fcanopy *\/ */
/* 			/\* Ls= Ls-Lssaved; /\* subtract the sunlit areas of all the upper layers saved in Lssaved; now we have sunlit area of current layer*\/ */
/* 			Lssaved = Ls+ Lssaved; /\* update the Lssaved for calculation of the next lower layer *\/ */
/* 			Ld=Ld-Ldsaved; /\* subtract the shaded areas of all the upper layers saved in Ldsaved; now we have shaded area of current layer*\/ */
/* 			Ldsaved=Ld+Ldsaved; /\* update the Ldsaved for calculation of the next lower layer *\/ */
/* 			Fsun=Ls/(Ls+Ld);   /\*fracction of sunlit leaves in the current layer *\/ */
/* 			Fshade=Ld/(Ls+Ld);  /\* fraction of shaded leaves in the current layer *\/ */
/*                         Itotal = Fsun * Isolar + Idiffuse; */
/* 		} */
		
/*                       /\* collecting the results *\/ */

/* 		layIdir[sp1++] = Isolar+Idiffuse;  /\* This is used to calculate sunlit leaf photosynthesis *\/ */
/* 		layIdiff[sp2++] = Idiffuse;       /\* ok, this is used to calculate shaded leaf photosynthesis *\/ */
/* 		layItotal[sp3++] = Itotal;        /\* this is used to evaluate evapotranspiration..Think about other factors that might influence energy balance to include *\/ */
/* 		layFsun[sp4++] = Fsun; */
/* 		layFshade[sp5++] = Fshade; */
/* 		layHeight[sp6++] = CumLAI / heightf; */
/* 	} */
/* } */


/* This is the function written by Joe Iverson */
void sunML(double Idir, double Idiff, double LAI, int nlayers,
	   double cosTheta, double kd, double chil, double heightf,
	   double maxIdir, double maxIdiff)
{
	extern int sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8;
	extern double layIdir[], layIdiff[], layItotal[], layFsun[], layFshade[], layHeight[], layMaxIdir[], layMaxIdiff[];
	double i;
	double k0, k1, k;
	double LAIi, CumLAI;
	double Isolar, Idiffuse, Ibeam, Iscat, Itotal, Iaverage, alphascatter;
	double maxIsolar, maxIdiffuse;
	double Ls, Ld;
	double Fsun, Fshade;
	alphascatter=0.8;
	k0 = sqrt(pow(chil ,2) + pow(tan(acos(cosTheta)),2));
	k1 = chil + 1.744*pow((chil+1.183),-0.733);
	k = k0/k1;
	if(k<0)
		k = -k;

	LAIi = LAI / nlayers;

	for(i=0;i<nlayers;i++)
	{
		CumLAI = LAIi * (i+0.5);
		
		Ibeam=Idir*cosTheta;
		Iscat = Ibeam * exp(-k *sqrt(alphascatter)* CumLAI)-Ibeam * exp(-k * CumLAI);
		
		
		Isolar = Ibeam * k;
		maxIsolar = maxIdir * cosTheta * k;
		Idiffuse = Idiff * exp(-kd * CumLAI) + Iscat;
		maxIdiffuse = maxIdiff * exp(-kd * CumLAI) + Iscat;
		
		Ls = (1-exp(-k*LAIi))*exp(-k*CumLAI)/k;
		Ld=LAIi-Ls;

		Fsun=Ls/(Ls+Ld);
		Fshade=Ld/(Ls+Ld);
		/*fraction intercepted*/
		Itotal =(Fsun*Isolar + Idiffuse) * (1-exp(-k*LAIi))/k;

		Iaverage =(Fsun*(Isolar + Idiffuse) + Fshade*Idiffuse) * (1-exp(-k*LAIi))/k;

		/* collecting the results */
		layIdir[sp1++] = Isolar + Idiffuse;
		layIdiff[sp2++] = Idiffuse;
		layItotal[sp3++] = Itotal;
		layMaxIdir[sp4++] = maxIsolar;
		layMaxIdiff[sp5++] = maxIdiffuse;
		layFsun[sp6++] = Fsun;
		layFshade[sp7++] = Fshade;
		layHeight[sp8++] = CumLAI/heightf;
	}
}



/* Additional Functions needed for EvapoTrans */

/* RH and Wind profile function */
void WINDprof(double WindSpeed, double LAI, int nlayers)
{
	int i;
	double k=0.7;
	double LI, CumLAI;
	double Wind;

	LI  = LAI / nlayers;
	for(i=0;i<nlayers;i++)
	{
		CumLAI = LI * (i + 1);
		Wind = WindSpeed * exp(-k * (CumLAI-LI));
		tmp3[tp3++] = Wind;
	}
}

void RHprof(double RH, int nlayers)
{
	int i;
	double kh, hsla, j;

	kh = 1 - RH;
	/* kh = 0.2; */
	/*kh = log(1/RH);*/
	for(i=0;i<nlayers;i++)
	{
		j = i + 1;
		hsla = RH * exp(kh * (j/nlayers));
		/*hsla = RH * exp(-kh * (j/nlayers));  new simpler version from Joe Iverson*/
		if(hsla > 1) hsla = 0.99; 
		tmp4[tp4++] = hsla;
	}
	/* It should return values in the 0-1 range */
}

void LNprof(double LeafN, double LAI, int nlayers, double kpLN)
{

	int i;
	double leafNla, LI, CumLAI;

	LI  = LAI / nlayers;
	for(i=0;i<nlayers;i++)
	{
		CumLAI = LI * (i + 1);
		leafNla = LeafN * exp(-kpLN * (CumLAI-LI));
		tmp5[tp5++] = leafNla;
	}

}

double TempToDdryA(double Temp)
{
	double DdryA;
	DdryA = 1.295163636 + -0.004258182 * Temp;
	return(DdryA); /* Units are kg m^-3 */
}

double TempToLHV(double Temp)
{
	double LHV;
	LHV = 2.501 + -0.002372727 * Temp;
	return(LHV); /* Units are MJ kg-1 */
}

double TempToSFS(double Temp)
{
	double SlopeFS;
	SlopeFS = 0.338376068 +  0.011435897 * Temp +  0.001111111 * pow(Temp,2);
	return(SlopeFS);
	/* The units are g m-3 K-1 */
}

double TempToSWVC(double Temp)
{
/* Temp should be in Celsius */
/* This is the arden buck equation */
	double SWVC;
	double a, b;
	a = (18.678 - Temp/234.5) * Temp;
	b = 257.14 + Temp;
	SWVC =  (6.1121 * exp(a/b))/10;
	return(SWVC); /* This is in hecto Pascals */
}

/* EvapoTrans function */
struct ET_Str EvapoTrans(double Rad, double Iave, double Imax, double Airtemperature, double RH,
			 double WindSpeed,double LeafAreaIndex, double CanopyHeight, double StomataWS, int ws,
			 double vmax2, double alpha2, double kparm, double theta, double beta, double Rd2, double b02, double b12)
{
	/* creating the structure to return */
	struct ET_Str tmp;
	struct c4_str tmpc4;

	const double LeafWidth = 0.04; /* This assumes a leaf 4 cm wide */
	const double kappa = 0.41;
	const double WindSpeedHeight = 2; /* This is the height at which the wind speed was measured */
	const double dCoef = 0.77;
	const double tau = 0.2;
	const double ZetaCoef = 0.026;
	const double ZetaMCoef = 0.13;
	const double LeafReflectance = 0.2;
	const double SpecificHeat = 1010; /* J kg-1 K-1 */
	const double StefanBoltzmann = 5.67037e-8; /* J m^-2 s^-1 K^-4 */

	double Tair, WindSpeedTopCanopy;
	double DdryA, LHV, SlopeFS, SWVC;
	double LayerRelativeHumidity, LayerWindSpeed, totalradiation;
	double LayerConductance, DeltaPVa, PsycParam, ga;
	double BoundaryLayerThickness, DiffCoef,LeafboundaryLayer;
	double d, Zeta, Zetam, ga0, ga1, ga2; 
	double ActualVaporPressure;
	double Ja, Ja2, Deltat;
	double PhiN, PhiN2;
	double TopValue, BottomValue;
	double EPen, TransR,EPries; 
	double OldDeltaT, rlc, ChangeInLeafTemp;
	double f; /* cloudiness function */
	double eprime; /* Apparent net clear sky emissivity */ 
	double rels;
	int Counter;

	WindSpeedTopCanopy = WindSpeed;
	Tair = Airtemperature;

	if(CanopyHeight < 0.1)
		CanopyHeight = 0.1; 

	DdryA = TempToDdryA(Tair);

	/* In the original code in WIMOVAC this is used in J kg-1
but Thornley and Johnson use it as MJ kg-1  */
	LHV = TempToLHV(Tair); /* This should be MJ kg-1*/
	LHV = LHV * 1e6; /* Now it is converted to Joules */
	SlopeFS = TempToSFS(Tair) * 1e-3;
	SWVC = TempToSWVC(Tair); /* this is hecto Pascals */
        /* Convert to kg/m3 */
	SWVC = (DdryA * 0.622 * SWVC)/1013.25; /* This last number is atmospheric pressure in hecto pascals */ 

	LayerRelativeHumidity = RH * 100;
	if(LayerRelativeHumidity > 100) 
		error("LayerRelativehumidity > 100"); 

        /* SOLAR RADIATION COMPONENT*/

        /* First calculate the Radiation term */
	/*' Convert light assuming 1 µmol PAR photons = 0.235 J */
        /* The next step converts from PAR photons to Joules */
        /* There are 2.35 x 10^5 Joules in a mol */
        /* or 0.235 Joules in a micro mol */
        /* A Watt is equivalent to J/s */
	totalradiation = Rad * 0.235; /* This is essentially Watts m^-2 */
        /* On a clear sky it may exceed 1000 in some parts of the world 
           Thornley and Johnson pg 400 */
/* This values can not possibly be higher than 650 */
        if(totalradiation > 650) error("total radiation too high");
        
	/* Ja = (2 * totalradiation * ((1 - LeafReflectance - tau) / (1 - tau))) * LeafAreaIndex; */
	/* It seems that it is not correct to multiply by the leaf area index. Thornley
	   and johnson pg. 400 suggest using (1-exp(-k*LAI) but note that for a full canopy
	   this is 1. so I'm using 1 as an approximation. */
	Ja = (2 * totalradiation * ((1 - LeafReflectance - tau) / (1 - tau)));

        /* The value below is only for leaf temperature */
	Ja2 = (2 * Iave * 0.235 * ((1 - LeafReflectance - tau) / (1 - tau)));

        /* Non iterative calculation of longwave radiation */
	ActualVaporPressure = (LayerRelativeHumidity/100) * SWVC;
	if(Imax < 1){
		rels = 0.3;
	}else{
		rels = Rad / (Imax * 0.8); /* This is an empirical way of detrmining how sunny it was */
	}
	if(rels > 1) rels = 1;
	if(rels < 0.3) rels = 0.3;

	f = 1.35 * rels - 0.35;
	if(f < 0.6) f = 0.6;

	eprime = (0.34 - 0.14 * sqrt(ActualVaporPressure));

	rlc = StefanBoltzmann * pow(273 + Airtemperature, 4) * eprime * f;

	/* Rprintf("rlc1 %.8f \n",rlc); */

        /* Net radiation */
	PhiN = Ja - rlc;

        /* PhiN 2 for leaf temperature */
	PhiN2 = Ja2 - rlc;

        /* AERODYNAMIC COMPONENT */
	Zeta = ZetaCoef * CanopyHeight;
	Zetam = ZetaMCoef * CanopyHeight;
	d = dCoef * CanopyHeight;

	if(WindSpeed < 0.5) WindSpeed = 0.5;

	LayerWindSpeed = WindSpeed;

	tmpc4 = c4photoC(Rad,Airtemperature,RH,vmax2,alpha2,kparm,theta,beta,Rd2,b02,b12,StomataWS, 380, ws); 

	LayerConductance = tmpc4.Gs;

        /* Convert from mmol H20/m2/s to m/s */
	LayerConductance = LayerConductance * (1.0/41000.0);
	/* LayerConductance = LayerConductance * 24.39 * 1e-6; Original from WIMOVAC */ 

	/* Thornley and Johnson use m s^-1 on page 418 */

	/* prevent errors due to extremely low Layer conductance */
	if(LayerConductance <=0)
		LayerConductance = 0.0001;

	if(SWVC < 0)
		error("SWVC < 0");
	/* Now RHprof returns relative humidity in the 0-1 range */
	DeltaPVa = SWVC * (1 - LayerRelativeHumidity / 100);

	/* Rprintf("DeltaPVa %.8f \n",DeltaPVa); */

	PsycParam =(DdryA * SpecificHeat) / LHV; /* This is in kg m-3 K-1 */

	/* Rprintf("PsycParam %.8f \n",PsycParam); */

	/* Calculation of ga */
	/* According to thornley and Johnson pg. 416 */
	ga0 = pow(kappa,2) * LayerWindSpeed;
	ga1 = log((WindSpeedHeight + Zeta - d)/Zeta);
	ga2 = log((WindSpeedHeight + Zetam - d)/Zetam);
	ga = ga0/(ga1*ga2);

	if(ga < 0)
		error("ga is less than zero");

	DiffCoef = (2.126 * 1e-5) + ((1.48 * 1e-7) * Airtemperature);
	BoundaryLayerThickness = 0.004 * sqrt(LeafWidth / LayerWindSpeed);
	LeafboundaryLayer = DiffCoef / BoundaryLayerThickness;

        /* In WIMOVAC this follows */
        /* ga = (ga * LeafboundaryLayer) / (ga + LeafboundaryLayer); */

	TopValue = PhiN2 * (1 / ga + 1 / LayerConductance) - LHV * DeltaPVa;
	BottomValue = LHV * (SlopeFS + PsycParam * (1 + ga / LayerConductance));
	Deltat = TopValue / BottomValue;


	/* This is the original from WIMOVAC*/
	Deltat = 0.01;
	ChangeInLeafTemp = 10;

	Counter = 0;
	while( (ChangeInLeafTemp > 0.5) && (Counter <= 10))
	{
		OldDeltaT = Deltat;

		rlc = 4 * (5.67*1e-8) * pow(273 + Airtemperature, 3) * Deltat;  

/* rlc=net long wave radiation emittted per second =radiation emitted per second - radiation absorbed per second=sigma*(Tair+deltaT)^4-sigma*Tair^4 */
 
/* Then you do a Taylor series about deltaT = 0 and keep only the zero and first order terms. */
 
/* or rlc=sigma*Tair^4+deltaT*(4*sigma*Tair^3)-sigma*Tair^4=4*sigma*Tair^3*deltaT */
 
/* where 4*sigma*Tair^3 is the derivative of sigma*(Tair+deltaT)^4 evaluated at deltaT=0, */

		PhiN2 = (Ja2 - rlc);

		TopValue = PhiN2 * (1 / ga + 1 / LayerConductance) - LHV * DeltaPVa;
		BottomValue = LHV * (SlopeFS + PsycParam * (1 + ga / LayerConductance));
		Deltat = TopValue / BottomValue;
		if(Deltat > 7)	Deltat = 7;
		if(Deltat < -7)	Deltat = -7;

		ChangeInLeafTemp = OldDeltaT - Deltat;
		if(ChangeInLeafTemp <0)
			ChangeInLeafTemp = -ChangeInLeafTemp;
		Counter++;
	}

        /* Net radiation */
	PhiN = Ja - rlc;

	if(PhiN < 0)
		PhiN = 0;

	/* Rprintf("Imax %.5f \n",Imax); */
	/* Rprintf("Deltat %.5f \n",Deltat); */
	/* Rprintf("rlc2 %.5f \n",rlc); */

	TransR = (SlopeFS * PhiN + (LHV * PsycParam * ga * DeltaPVa)) / (LHV * (SlopeFS + PsycParam * (1 + ga / LayerConductance)));

	EPen = (((SlopeFS * PhiN) + LHV * PsycParam * ga * DeltaPVa)) / (LHV * (SlopeFS + PsycParam));

	EPries = 1.26 * ((SlopeFS * PhiN) / (LHV * (SlopeFS + PsycParam)));

	/* Rprintf("PhiN %.8f \n", PhiN); */
	/* Rprintf("Ja %.8f \n", Ja); */
	/* Rprintf("SlopeFS %.8f \n", SlopeFS); */
	/* Rprintf("LHV %.8f \n", LHV); */
	/* Rprintf("PsycParam %.8f \n", PsycParam); */
	/* Rprintf("DeltaPVa %.8f \n", DeltaPVa); */
	/* Rprintf("ga %.8f \n", ga); */
	/* Rprintf("LayerConductance %.8f \n", LayerConductance); */

        /* For now */
	TransR = EPen;

	/* This values need to be converted from Kg/m2/s to
	   mmol H20 /m2/s according to S Humphries */
	/* 1e3 - kgrams to grams  */
	/* 1e3 - mols to mmols */
        /* grams to mols - 18g in a mol */
	/* Let us return the structure now */

	tmp.TransR = TransR * 1e6 / 18; 
	tmp.EPenman = EPen * 1e6 / 18; 
	tmp.EPriestly = EPries * 1e6 / 18; 
	tmp.Deltat = Deltat;
	tmp.LayerCond = LayerConductance * 41000;   
	/*    tmp.LayerCond = RH2;   */
	/*   tmp.LayerCond = 0.7; */
	return(tmp);
}

/* Soil Evaporation Function */
/* Variables I need */
/* LAI = Leaf Area Index */
/* k = extinction coefficient */
/* Temp = Air Temperature */
/* DirectRad = Direct Total Radiation */
/* awc, wiltp, fieldc = available water content, wilting point and field capacity */
/* winds = wind speed */

double SoilEvapo(double LAI, double k, double AirTemp, double IRad,
		 double awc, double fieldc, double wiltp, double winds, double RelH, double rsec,
                 int soiltype){

	struct soilText_str soTexS;

	double SoilArea;
	double SoilTemp;
	double Up; /*Maximum Dimensionless Uptake Rate */
	double TotalRadiation;
	double BoundaryLayerThickness, DiffCoef;
	double SoilBoundaryLayer, Ja, rlc;
	double PhiN, PsycParam, DeltaPVa;
	double Evaporation = 0.0;  
	double DdryA, LHV, SlopeFS, SWVC;

	double rawc; /* Relative available water content */

	int method = 1;

	/* some constants */
	const double SoilClodSize = 0.04;
	const double SoilReflectance = 0.2;
	const double SoilTransmission = 0.01;
	const double SpecificHeat = 1010;
	const double StefanBoltzman = 5.67e-8; /* Stefan Boltzman Constant */

	const double cf2 = 3600 * 1e-3 * 18 * 1e-6 * 10000; 
        /* Do I understand this conversion factor ?*/ 
        /* 3600 converts seconds to hour */

	/* Specify the soil type */
	soTexS = soilTchoose(soiltype);

	if(fieldc < 0){
		fieldc = soTexS.fieldc;
	}
	if(wiltp < 0){
		wiltp = soTexS.wiltp;
	}

	/* For Transpiration */
	/* 3600 converts seconds to hours */
	/* 1e-3 converts mili mols to mols */
	/* 18 is the grams in one mol of H20 */
	/* 1e-6 converts g to Mg */
	/* 10000 scales from meter squared to hectare */

	/* Let us assume a simple way of calculating the proportion of the
	   soil with direct radiation */
	SoilArea = exp(-k * LAI);

	/* For now the temperature of the soil will be the same as the air.
	   At a later time this can be made more accurate. I looked at the
	   equations for this and the issue is that it is strongly dependent on
	   depth. Since the soil model now has a single layer, this cannot be
	   implemented correctly at the moment.  */

	SoilTemp = AirTemp;

	rawc = (awc - wiltp) / (fieldc - wiltp);
        if(rawc < 1e-7) rawc = 1e-7;
	if(rawc > 1) rawc = 1;

	/* Up adjusts for soil drying */
	if(awc > fieldc){
		Up =  1;
	}else{
		Up = exp(rawc * 5) / exp(5);
	}  
	/* This value will be close to 1 with saturated soil and close to zero as the soil dries*/
        /* This is an empirical relationship that I made up */

	/* Total Radiation */
	/*' Convert light assuming 1 µmol PAR photons = 0.235 J/s Watts*/
	/* At the moment soil evaporation is grossly overestimated. In WIMOVAC
	   the light reaching the last layer of leaves is used. Here instead
	   of calculating this again, I will for now assume a 10% as a rough
	   estimate. Note that I could maybe get this since layIdir and
	   layIDiff in sunML are external variables.  Rprintf("IRad
	   %.5f",layIdir[0],"\n"); Update: 03-13-2009. I tried printing this
	   value but it is still too high and will likely overestimate soil
	   evaporation. However, this is still a work in progress.
	*/
	IRad *= rsec;  /* Radiation soil evaporation coefficient  */

	TotalRadiation = IRad * 0.235;
 
	DdryA = TempToDdryA(AirTemp);
	LHV = TempToLHV(AirTemp) * 1e6 ; 
/* Here LHV is given in MJ kg-1 and this needs to be converted
   to Joules kg-1  */
	SlopeFS = TempToSFS(AirTemp) * 1e-3; /* Need to change from g/m3/K to kg/m3/K */
	SWVC = TempToSWVC(AirTemp) * 1e-3; /* Need to change from g/m3 to kg/m3 */

	PsycParam = (DdryA * SpecificHeat) / LHV;
	DeltaPVa = SWVC * (1 - RelH);

	BoundaryLayerThickness = 4e-3 * sqrt(SoilClodSize / winds); 
	DiffCoef = 2.126e-5 * 1.48e-7 * SoilTemp;
	SoilBoundaryLayer = DiffCoef / BoundaryLayerThickness;

	Ja = 2 * TotalRadiation * ((1 - SoilReflectance - SoilTransmission) / (1 - SoilTransmission));

	rlc = 4 * StefanBoltzman * pow((273 + SoilTemp),3) * 0.005;
/* the last term should be the difference between air temperature and
   soil. This is not actually calculated at the moment. Since this is
   mostly relevant to the first soil layer where the temperatures are
   similar. I will leave it like this for now. */

	PhiN = Ja - rlc; /* Calculate the net radiation balance*/
	if(PhiN < 0) PhiN = 1e-7;

	/* Rprintf("PhiN: %.6f \n", PhiN);  */
	/* Rprintf("Ja: %.6f \n", Ja);  */
	/* Rprintf("rlc: %.6f \n", rlc);  */
	/* Rprintf("total radiation: %.6f \n", TotalRadiation);  */
	/* Rprintf("IRad: %.6f \n", IRad);  */
	/* Rprintf("SlopeFS: %.6f \n", SlopeFS); */
	/* Rprintf("LHV: %.6f \n", LHV);   */
	/* Rprintf("PsycParam: %.6f \n", PsycParam);  */

	/* Priestly-Taylor */
	if(method == 0){
		Evaporation = 1.26 * (SlopeFS * PhiN) / (LHV * (SlopeFS + PsycParam));
	}else{
		/* Penman-Monteith */
		Evaporation = (SlopeFS * PhiN + LHV * PsycParam * SoilBoundaryLayer * DeltaPVa) / (LHV * (SlopeFS + PsycParam));
	}

/* What are the actual units of Evaporation? Below I assume these are kg H20/m2/s*/


/*  Report back the soil evaporation rate in Units mmoles/m2/s */
/*     Evaporation = Evaporation * 1000:   ' Convert Kg H20/m2/s to g H20/m2/s */
/*     Evaporation = Evaporation / 18:     ' Convert g H20/m2/s to moles H20/m2/s */
/*     Evaporation = Evaporation * 1000:   ' Convert moles H20/m2/s to mmoles H20/m2/s */

/* Assuming that it returns data on kg H20/m2/s I will convert to mm/hr */
        /* Evaporation *= 1e-3 * 3600; */

	Evaporation *= 1e6/18;

	/* Adding the area dependence and the effect of drying */
	/* Converting from m2 to ha (times 1e4) */
	/* Converting to hour */
	Evaporation *= SoilArea * Up * cf2; 

	if(Evaporation < 0) Evaporation = 1e-6;

	return(Evaporation); /* The units returned are Mg/ha/hr */
}


/* CanA function */
/* This file will contain the function which will consolidate
   the previous functions into one. The idea is to go from
   existing LAI and weather conditions to canopy assimilation
   and transpiration */

struct Can_Str CanAC(double LAI,int DOY, int hr,double solarR,double Temp,
	             double RH,double WindSpeed,double lat,int nlayers, double Vmax,
		     double Alpha, double Kparm, double theta, double beta,
		     double Rd, double Catm, double b0, double b1,
                     double StomataWS, int ws, double kd, double chil, double heightf,
		     double leafN, double kpLN, double lnb0, double lnb1, int lnfun)
{

	struct ET_Str tmp5_ET, tmp6_ET;
	struct Can_Str ans;
	struct c4_str tmpc4;
	struct c4_str tmpc42;

	int i;
	double Idir, Idiff, cosTh, maxIdir, maxIdiff;
	double maxIDir, maxIDiff;
	double LAIc;
	double IDir, IDiff, Itot, rh, WS;
	double pLeafsun, pLeafshade;
	double Leafsun, Leafshade;
	double CanHeight;

	double vmax1, leafN_lay;
	double TempIdir,TempIdiff,AssIdir,AssIdiff;

	double CanopyA, CanopyT;

	const double cf = 3600 * 1e-6 * 30 * 1e-6 * 10000;
	const double cf2 = 3600 * 1e-3 * 18 * 1e-6 * 10000; 

	/* For Assimilation */
	/* 3600 converts seconds to hours */
	/* 1e-6 converts micro mols to mols */
	/* 30 is the grams in one mol of CO2 */
	/* 1e-6 converts g to Mg */
	/* 10000 scales from meter squared to hectare */

	/* For Transpiration */
	/* 3600 converts seconds to hours */
	/* 1e-3 converts mili mols to mols */
	/* 18 is the grams in one mol of H20 */
	/* 1e-6 converts g to Mg */
	/* 10000 scales from meter squared to hectare */

	lightME(lat,DOY,hr);

	Idir = tmp1[0] * solarR;
	Idiff = tmp1[1] * solarR;
	cosTh = tmp1[2];
	maxIdir = tmp1[3];
	maxIdiff = tmp1[4];
    
	sunML(Idir,Idiff,LAI,nlayers,cosTh, kd, chil, heightf, maxIdir, maxIdiff);

	/* results from multilayer model */
	LAIc = LAI / nlayers;
	/* Next I need the RH and wind profile */
	RHprof(RH,nlayers);
	WINDprof(WindSpeed,LAI,nlayers);

	LNprof(leafN, LAI, nlayers, kpLN);
	/* It populates tmp5 */

	/* Next use the EvapoTrans function */
	CanopyA=0.0;
	CanopyT=0.0;
	for(i=0;i<nlayers;i++)
	{
		leafN_lay = tmp5[--tp5];
		if(lnfun == 0){
			vmax1 = Vmax;
		}else{
			vmax1 = leafN_lay * lnb1 + lnb0;
               /* For now alpha is not affected by leaf nitrogen */
		}

		IDir = layIdir[--sp1];
		Itot = layItotal[--sp3];

		maxIDir = layMaxIdir[--sp4];
		maxIDiff = layMaxIdiff[--sp5];

		rh = tmp4[--tp4];
		WS = tmp3[--tp3];
		pLeafsun = layFsun[--sp6];
		CanHeight = layHeight[--sp8];
		Leafsun = LAIc * pLeafsun;
		tmp5_ET = EvapoTrans(IDir,Itot,maxIDir,Temp,rh,WS,LAIc,CanHeight,StomataWS,ws,vmax1,Alpha,Kparm,theta,beta,Rd,b0,b1);
		TempIdir = Temp + tmp5_ET.Deltat;
		tmpc4 = c4photoC(IDir,TempIdir,rh,vmax1,Alpha,Kparm,theta,beta,Rd,b0,b1,StomataWS, Catm, ws);
		AssIdir = tmpc4.Assim;


		IDiff = layIdiff[--sp2];
		pLeafshade = layFshade[--sp7];
		Leafshade = LAIc * pLeafshade;
		tmp6_ET = EvapoTrans(IDiff,Itot,maxIDiff,Temp,rh,WS,LAIc,CanHeight,StomataWS,ws,vmax1,Alpha,Kparm,theta,beta,Rd,b0,b1);
		TempIdiff = Temp + tmp6_ET.Deltat;
		tmpc42 = c4photoC(IDiff,TempIdiff,rh,vmax1,Alpha,Kparm,theta,beta,Rd,b0,b1,StomataWS, Catm, ws);
		AssIdiff = tmpc42.Assim;

		CanopyA += Leafsun * AssIdir + Leafshade * AssIdiff;
		CanopyT += Leafsun * tmp5_ET.TransR + Leafshade * tmp6_ET.TransR;
	}
	/*## These are micro mols of CO2 per m2 per sec for Assimilation
	  ## and mili mols of H2O per m2 per sec for Transpiration
	  ## Need to convert to 
	  ## 3600 converts seconds to hours
	  ## 10^-6 converts micro mols to mols
	  ## 30 converts mols of CO2 to grams
	  ## (1/10^6) converts grams to Mg
	  ## 10000 scales up to ha */
/* A similar conversion is made for water but
   replacing 30 by 18 and mili mols are converted to
   mols (instead of micro) */
	ans.Assim = cf * CanopyA ;
/* CanopyTrans can apparently go crazy */
	/* Rprintf("CanopyT %.8f \n",CanopyT); */
	if(CanopyT > 20){
		Rprintf("CanopyT %.5f \n", CanopyT);
		error("CanopyT is too high");
	}
	ans.Trans = cf2 * CanopyT; 
	return(ans);
}


/* This is a new function that attempts to keep a water budget and then
   calcualte an empirical coefficient that reduces the specific leaf area.
   This results from the general idea that water stress reduces first the 
   rate of leaf expansion. */ 

/* This is meant to be a simple function that calculates a
   simple empirical coefficient that reduces specifi leaf area
   according to the water stress of the plant. This is done
   for now, with a very simple empirical approach. */

struct ws_str watstr(double precipit, double evapo, double cws, double soildepth, 
                     double fieldc, double wiltp, double phi1, double phi2, 
		     double smthresh,
		     int soiltype, /* soil type indicator */ 
		     int wsFun) /* flag for which water stress function to use */
{

	struct ws_str tmp;
	struct soilText_str soTexS;
	const double g = 9.8; /* m / s-2  ##  http://en.wikipedia.org/wiki/Standard_gravity */
	/* Variables */
	double precipM;
	/* available water and per hectare */
	double aw, naw, raw, theta_s; 
	double K_psim, J_w;
	double pawha, Newpawha, npaw;
	double runoff = 0.0, runoff2 = 0.0, drainage = 0.0;
	/* variable needed for calculation of water stress*/
	double wsPhoto = 0.0, wsSpleaf, phi10;
	double slp = 0.0, intcpt = 0.0, theta = 0.0; 
	double Nleach = 0.0;
	/* Nleach is the NO3 leached and Ts is the sand content of the soil*/

	/* Specify the soil type */
	soTexS = soilTchoose(soiltype);
/*   Ts = soTexS.sand; */

	if(fieldc < 0){
		fieldc = soTexS.fieldc;
	}
	if(wiltp < 0){
		wiltp = soTexS.wiltp;
	}

	theta_s = soTexS.satur;

	/* unit conversion for precip */
	precipM = precipit * 1e-3; /* convert precip in mm to m*/

	/*    cws is current water status */
	/*    available water */

	aw = precipM + cws;

/* These equations are not correct as runoff would only occur when it exceeds
   saturation, but from the point of view of a crop only field capacity matters */
/* I'm not sure about what to do about this */

	if(aw > theta_s){ 
		runoff = aw - theta_s; /* Here runoff is interpreted as water content exceeding saturation level */
		/* Need to convert to units used in the Parton et al 1988 paper. */
		/* The data comes in mm/hr and it needs to be in cm/month */
		runoff2 = runoff * 0.10 * (1/24*30);
		Nleach = runoff /18 * (0.2 + 0.7 * soTexS.sand);
		aw = theta_s;
	}


	/* Tipping bucket need to collect it if want to estimate runoff */ 
	/* plant available water per hectare (pawha) */
	pawha = (aw - wiltp) * 1e4 * soildepth;
	/* The density of water is 998.2 kg/m3 at 20 degrees Celsius */
	/* or 0.9882 Mg/m3 */
	/* pawha is plant available water (m3) per hectare */
	/* evapo is demanded water (Mg) per hectare */

	Newpawha = pawha - evapo / 0.9882; /* New version 04-27-2009 */

	/*  Here both are in m3 of water per ha-1 so this */
	/*  subtraction should be correct */
	/* go back to original units of water in the profile */

	npaw = Newpawha * 1e-4 * (1/soildepth); /* New 04-27-2009 */

	if(npaw < 0) npaw = 0.0;

	naw = npaw + wiltp;

        /* Calculating the soil water potential based on equations from Norman and Campbell */
	/* tmp.psim = soTexS.air_entry * pow((naw/soTexS.fieldc*1.1),-soTexS.b) ; */
	/* New version of the soil water potential is based on
	 * "Dynamic Simulation of Water Deficit Effects upon Maize
	 * Yield" R. F. Grant Agricultural Systems. 33(1990) 13-39. */
        tmp.psim = -exp(log(0.033) + ((log(fieldc) - log(naw))/(log(fieldc) - log(wiltp)) * (log(1.5) - log(0.033)))) * 1e3; /* This last term converts from MPa to kPa */

	/* This is drainage */
	if(naw > fieldc){
	  K_psim = soTexS.Ks * pow((soTexS.air_entry/tmp.psim),2+3/soTexS.b); /* This is hydraulic conductivity */
	  J_w = -K_psim * (-tmp.psim/(soildepth*0.25)) - g * K_psim ; /*  Campbell, pg 129 do not ignore the graviational effect. I multiply soil depth by 0.5 to calculate the average depth*/
	  drainage = J_w * 3600 * 0.9882 * 1e-3; /* This is flow in m3 / (m^2 * hr). */
	  naw = naw + drainage;
	}

	/* three different type of equations for modeling the effect of water stress on vmax and leaf area expansion. 
	   The equation for leaf area expansion is more severe than the one for vmax. */
	if(wsFun == 0){ /* linear */
		slp = 1/(fieldc - wiltp);
		intcpt = 1 - fieldc * slp;
		wsPhoto = slp * naw + intcpt ;
	}else
	if(wsFun == 1){
		phi10 = (fieldc + wiltp)/2;
		wsPhoto = 1/(1 + exp((phi10 - naw)/ phi1));
	}else
	if(wsFun == 2){
		slp = (1 - wiltp)/(fieldc - wiltp);
		intcpt = 1 - fieldc * slp;
		theta = slp * naw + intcpt ;
		wsPhoto = (1 - exp(-2.5 * (theta - wiltp)/(1 - wiltp))) / (1 - exp(-2.5));
	}else
	if(wsFun == 3){
		wsPhoto = 1;
	}else
	if(wsFun == 5){
		raw = (naw - wiltp)/(fieldc - wiltp);
		if(raw > smthresh){
			wsPhoto = 1;
		}else{
			wsPhoto = raw / smthresh;
		}
	}


	if(wsPhoto <= 0 )
		wsPhoto = 1e-20; /* This can be mathematically lower than zero in some cases but I should prevent that. */

	wsSpleaf = pow(naw,phi2) * 1/pow(fieldc,phi2); 
	if(wsFun == 3){ 
		wsSpleaf = 1;
	}

/* Apparently wsPhoto and wsSpleaf can be greater than 1 */
        if(wsPhoto > 1) wsPhoto = 1;
	if(wsSpleaf > 1) wsSpleaf = 1;

	/* returning the structure*/
	tmp.rcoefPhoto = wsPhoto;
	tmp.rcoefSpleaf = wsSpleaf;
	tmp.awc = naw;
	tmp.runoff = runoff;
	tmp.Nleach = Nleach;
	tmp.drainage = drainage;
	return(tmp);
}

/* Function to simulate the multilayer behavior of soil water. In the
   future this could be coupled with Campbell (BASIC) ideas to
   esitmate water potential. */
struct soilML_str soilML(double precipit, double transp, double *cws, double soildepth, double *depths, double fieldc, double wiltp, double phi1, double phi2, double smthresh, struct soilText_str soTexS, int wsFun, int layers, double rootDB, double LAI, double k, double AirTemp, double IRad, double winds, double RelH, int hydrDist, double rfl, double rsec, double rsdf){

	struct rd_str tmp4;
	struct seqRD_str tmp3;
	struct soilML_str tmp;
        /* Constant */
	/* const double G = 6.67428e-11;  m3 / (kg * s-2)  ##  http://en.wikipedia.org/wiki/Gravitational_constant */
	const double g = 9.8; /* m / s-2  ##  http://en.wikipedia.org/wiki/Standard_gravity */
	/* Variables */
	double waterIn, oldWaterIn = 0.0;
/* Here is a convention aw is available water in volume and awc
   is available water content as a fraction of the soil section being investigated.
   paw is plant available water aw - wiltp */
	double aw, paw, awc, awc2, Newpawha, raw;
	double drainage = 0.0;
	double wsPhoto = 0.0, wsSpleaf = 0.0, phi10;
	double wsPhotoCol = 0.0, wsSpleafCol = 0.0;
	double slp = 0.0, intcpt = 0.0, theta = 0.0; 
	double Nleach = 0.0;
	double layerDepth;
	double diffw;
	double rootATdepth, rootDepth;
	double EvapoTra = 0.0, oldEvapoTra = 0.0, Sevap = 0.0, Ctransp = 0.0;
	double psim1 = 0.0, psim2 = 0.0, K_psim = 0.0, J_w = 0.0, dPsim = 0.0;
	double theta_s; /* This is the saturated soil water content. Larger than FieldC.*/
	int i;
	int j = layers - 1; 

	/* Specify the soil type */

	if(fieldc < 0){
		fieldc = soTexS.fieldc;
	}
	if(wiltp < 0){
		wiltp = soTexS.wiltp;
	}

	theta_s = soTexS.satur;
	/* rooting depth */
	/* Crude empirical relationship between root biomass and rooting depth*/
	rootDepth = rootDB * rsdf;
	if(rootDepth > soildepth) rootDepth = soildepth;

	tmp3 = seqRootDepth(rootDepth,layers);
	tmp4 = rootDist(layers,rootDepth,&depths[0],rfl);

	/* unit conversion for precip */
	waterIn = precipit * 1e-3; /* convert precip in mm to m*/

	for(j=0,i=layers-1;j<layers;j++,i--){
	/* for(i=0;i<layers;i++){ */
		/* It decreases because I increase the water content due to precipitation in the last layer first*/

		/* This supports unequal depths. */
		if(i == 0){
			layerDepth = depths[1];
		}else{
			layerDepth = depths[i] - depths[i-1];
		}


		if(hydrDist > 0){
			/* For this section see Campbell and Norman "Environmental BioPhysics" Chapter 9*/
			/* First compute the matric potential */
			psim1 = soTexS.air_entry * pow((cws[i]/theta_s),-soTexS.b) ; /* This is matric potential of current layer */
			if(i > 0){
				psim2 = soTexS.air_entry * pow((cws[i-1]/theta_s),-soTexS.b) ; /* This is matric potential of next layer */
				dPsim = psim1 - psim2;
				/* The substraction is from the layer i - (i-1). If this last term is positive then it will move upwards. If it is negative it will move downwards. Presumably this term is almost always positive. */
			}else{
				dPsim = 0;
			}
			K_psim = soTexS.Ks * pow((soTexS.air_entry/psim1),2+3/soTexS.b); /* This is hydraulic conductivity */
			J_w = K_psim * (dPsim/layerDepth) - g * K_psim ; /*  Campbell, pg 129 do not ignore the graviational effect*/
                        /* Notice that K_psim is positive because my
                            reference system is reversed */
			/* This last result should be in kg/(m2 * s)*/
			 J_w *= 3600 * 0.9882 * 1e-3 ; /* This is flow in m3 / (m^2 * hr). */
			/* Rprintf("J_w %.10f \n",J_w);  */
			if(i == (layers-1) && J_w < 0){
					/* cws[i] = cws[i] + J_w /
					 * layerDepth; Although this
					 * should be done it drains
					 * the last layer too much.*/
					drainage += J_w;
			}else{
				if(i > 0){
					cws[i] = cws[i] -  J_w / layerDepth;
					cws[i - 1] =  cws[i-1] +  J_w / layerDepth;
				}else{
					cws[i] = cws[i] -  J_w / layerDepth;
				}
			}
		}

		 if(cws[i] > theta_s) cws[i] = theta_s; 
		/* if(cws[i+1] > fieldc) cws[i+1] = fieldc; */
		 if(cws[i] < wiltp) cws[i] = wiltp; 
		/* if(cws[i+1] < wiltp) cws[i+1] = wiltp;  */

		aw = cws[i] * layerDepth;
/* Available water (for this layer) is the current water status times the layer depth */

		if(waterIn > 0){
			/* There is some rain. Need to add it.*/
			aw += waterIn / layers + oldWaterIn; /* They are both in meters so it works */
                        /* Adding the same amount to water to each layer */
                        /* In case there is overflow */
			/* diffw = fieldc * layerDepth - aw; */
			diffw = theta_s * layerDepth - aw;

			if(diffw < 0){
				/* This means that precipitation exceeded the capacity of the first layer */
				/* Save this amount of water for the next layer */
				oldWaterIn = -diffw;
				/* aw = fieldc * layerDepth; */
				aw = theta_s * layerDepth;
			}else{
				oldWaterIn = 0.0;
			}
		}

		/* Root Biomass */
		rootATdepth = rootDB * tmp4.rootDist[i];
		tmp.rootDist[i] = rootATdepth;
/* Plant available water is only between current water status and permanent wilting point */
		/* Plant available water */
		paw = aw - wiltp * layerDepth;
		if(paw < 0) paw = 0; 

		if(i == 0){
			/* Only the first layer is affected by soil evaporation */
			awc2 = aw / layerDepth;
			/* SoilEvapo function needs soil water content  */
			Sevap = SoilEvapo(LAI,k,AirTemp,IRad,awc2,fieldc,wiltp,winds,RelH,rsec, 6);
                        /* Important note the 6 at the end is arbitrary and it won't affect the calculations */
			/* I assume that crop transpiration is distributed simlarly to
			   root density.  In other words the crop takes up water proportionally
			   to the amount of root in each respective layer.*/
			Ctransp = transp*tmp4.rootDist[0];
			EvapoTra = Ctransp + Sevap;
			Newpawha = (paw * 1e4) - EvapoTra / 0.9982; /* See the watstr function for this last number 0.9882 */
			/* The first term in the rhs (paw * 1e4) is the m3 of water available in this layer.
			   EvapoTra is the Mg H2O ha-1 of transpired and evaporated water. 1/0.9882 converts from Mg to m3 */
		}else{
			Ctransp = transp*tmp4.rootDist[i];
			EvapoTra = Ctransp;
			Newpawha = (paw * 1e4) - (EvapoTra + oldEvapoTra);
		}

		if(Newpawha < 0){
/* If the Demand is not satisfied by this layer. This will be stored and added to subsequent layers*/
			oldEvapoTra = -Newpawha;
			 aw = wiltp * layerDepth; 
		}

		paw = Newpawha / 1e4 ;
		awc = paw / layerDepth + wiltp;   

/* This might look like a weird place to populate the structure, but is more convenient*/
		tmp.cws[i] = awc;

		if(wsFun == 0){
			slp = 1/(fieldc - wiltp);
			intcpt = 1 - fieldc * slp;
			wsPhoto = slp * awc + intcpt ;
		}else
		if(wsFun == 1){
			phi10 = (fieldc + wiltp)/2;
			wsPhoto = 1/(1 + exp((phi10 - awc)/ phi1));
		}else
		if(wsFun == 2){
			slp = (1 - wiltp)/(fieldc - wiltp);
			intcpt = 1 - fieldc * slp;
			theta = slp * awc + intcpt ;
			wsPhoto = (1 - exp(-2.5 * (theta - wiltp)/(1 - wiltp))) / (1 - exp(-2.5));
		}else
		if(wsFun == 3){
			wsPhoto = 1;
		}else
		if(wsFun == 5){
			raw = (awc - wiltp)/(fieldc - wiltp);
		        if(raw > smthresh){
			   wsPhoto = 1;
		        }else{
			   wsPhoto = raw / smthresh;
			}
		}



		if(wsPhoto <= 0 )
			wsPhoto = 1e-20; /* This can be mathematically lower than zero in some cases but I should prevent that. */

		wsPhotoCol += wsPhoto;

		wsSpleaf = pow(awc,phi2) * 1/pow(fieldc,phi2); 
		if(wsFun == 3){ 
			wsSpleaf = 1;
		}
		wsSpleafCol += wsSpleaf;

	}

	if(waterIn > 0){ 
		drainage = waterIn;
		/* Need to convert to units used in the Parton et al 1988 paper. */
		/* The data comes in mm/hr and it needs to be in cm/month */
		Nleach = drainage * 0.1 * (1/24*30) / (18 * (0.2 + 0.7 * soTexS.sand));
	}

/* Apparently wsPhoto and wsSpleaf can be greater than 1 */
        if(wsPhoto > 1) wsPhoto = 1;
	if(wsSpleaf > 1) wsSpleaf = 1;

/* returning the structure */
	tmp.rcoefPhoto = (wsPhotoCol/layers);
	tmp.drainage = drainage;
	tmp.Nleach = Nleach;
	tmp.rcoefSpleaf = (wsSpleafCol/layers);
	tmp.SoilEvapo = Sevap;

	return(tmp);
}

/* Respiration. It is assumed that some of the energy produced by the
   plant has to be used in basic tissue maintenance. This is handled
   quite empirically by some relationships developed by McCree (1970)
   and Penning de Vries (1972) */

double resp(double comp, double mrc, double temp){

	double ans;

	ans = comp *  (1 - (mrc * pow(2,(temp/10.0))));

	if(ans <0) ans = 0;

	return(ans);

}

/* Function to select the correct dry biomass partitioning coefficients */
/* It should take a vector of size 24 as an argument and return a structure with four numbers */

struct dbp_str sel_dbp_coef(double coefs[25], double TherPrds[6], double TherTime){

	struct dbp_str tmp;

	tmp.kLeaf = 0.0;
	tmp.kStem = 0.0;
	tmp.kRoot = 0.0;
	tmp.kRhiz = 0.0;
	tmp.kGrain = 0.0; /* kGrain is always zero except for the last thermal period */


	if(TherTime < TherPrds[0])
	{
		tmp.kStem = coefs[0];
		tmp.kLeaf = coefs[1];
		tmp.kRoot = coefs[2];
		tmp.kRhiz = coefs[3];

	} else
		if( TherTime < TherPrds[1] )
		{
			tmp.kStem = coefs[4];
			tmp.kLeaf = coefs[5];
			tmp.kRoot = coefs[6];
			tmp.kRhiz = coefs[7];

		} else
			if( TherTime < TherPrds[2])
			{
				tmp.kStem = coefs[8];
				tmp.kLeaf = coefs[9];
				tmp.kRoot = coefs[10];
				tmp.kRhiz = coefs[11];

			} else
				if(TherTime < TherPrds[3])
				{
					tmp.kStem = coefs[12];
					tmp.kLeaf = coefs[13];
					tmp.kRoot = coefs[14];
					tmp.kRhiz = coefs[15];

				} else
					if(TherTime < TherPrds[4])
					{
						tmp.kStem = coefs[16];
						tmp.kLeaf = coefs[17];
						tmp.kRoot = coefs[18];
						tmp.kRhiz = coefs[19];

					} else
						if(TherTime < TherPrds[5])
						{
							tmp.kStem = coefs[20];
							tmp.kLeaf = coefs[21];
							tmp.kRoot = coefs[22];
							tmp.kRhiz = coefs[23];
							tmp.kGrain = coefs[24];
						}else{
							Rprintf("TherPrds[5]: %.1f TherTime %.1f \n",TherPrds[5],TherTime);
							error("Thermal time larger than thermal periods");
						}


	return(tmp);

}

struct seqRD_str seqRootDepth(double to, int lengthOut){

	struct seqRD_str tmp;
	int i;
	double by;

	/* This is because in this case from is always zero */
	by = to / lengthOut;

	for(i=0;i<=lengthOut;i++){

		tmp.rootDepths[i] = i * by;

	}

	return(tmp);

}


struct rd_str rootDist(int layer, double rootDepth, double *depthsp, double rfl){

	struct rd_str tmp;  
	int i, j, k;
	double layerDepth = 0.0;
	double CumLayerDepth = 0.0;
	double CumRootDist = 1.0;
	double rootDist[layer];
	double ca = 0.0, a = 0.0;

	for(i=0;i<layer;i++){

		if(i == 0){
			layerDepth = depthsp[1];
		}else{
			layerDepth = depthsp[i] - depthsp[i-1];
		}

		CumLayerDepth += layerDepth;

		if(rootDepth > CumLayerDepth){
			CumRootDist++;
		}
	}

	for(j=0;j<layer;j++){
		if(j < CumRootDist){ 
			a = dpois(j+1,CumRootDist*rfl,0);
			rootDist[j] = a;
			ca += a;
		}else{
			rootDist[j] = 0;
		}
	}

	for(k=0;k<layer;k++){
		tmp.rootDist[k] = rootDist[k] / ca; 
	}

	return(tmp);
}


struct soilText_str soilTchoose(int soiltype){

	/* This function is based on Campbell and Norman.
	   Introduction to Environmental Biophysics. pg 130. */

	struct soilText_str tmp;

	tmp.silt = 0;
	tmp.clay = 0;
	tmp.sand = 0;
	tmp.air_entry = 0;
	tmp.b = 0;
	tmp.Ks = 0;
	tmp.satur = 0;
	tmp.fieldc = 0;
	tmp.wiltp = 0;

	if(soiltype == 0){
	/* sand soil */
	tmp.silt = 0.05;
	tmp.clay = 0.03;
	tmp.sand = 0.92;
	tmp.air_entry = -0.7;
	tmp.b = 1.7;
	tmp.Ks = 5.8e-3;
	tmp.satur = 0.87;
	tmp.fieldc = 0.09;
	tmp.wiltp = 0.03;

	} else

	if(soiltype == 1){
	/* loamy sand */
	tmp.silt = 0.12;
	tmp.clay = 0.07;
	tmp.sand = 0.81;
	tmp.air_entry = -0.9;
	tmp.b = 2.1;
	tmp.Ks = 1.7e-3;
	tmp.satur = 0.72;
	tmp.fieldc = 0.13;
	tmp.wiltp = 0.06;

	} else

	if(soiltype == 2){
	/* sandy loam */
	tmp.silt = 0.25;
	tmp.clay = 0.10;
	tmp.sand = 0.65;
	tmp.air_entry = -1.5;
	tmp.b = 3.1;
	tmp.Ks = 7.2e-4;
	tmp.satur = 0.57;
	tmp.fieldc = 0.21;
	tmp.wiltp = 0.10;

	} else

	if(soiltype == 3){
	/* loam */
	tmp.silt = 0.40;
	tmp.clay = 0.18;
	tmp.sand = 0.52;
	tmp.air_entry = -1.1;
	tmp.b = 4.5;
	tmp.Ks = 3.7e-4;
	tmp.satur = 0.57;
	tmp.fieldc = 0.27;
	tmp.wiltp = 0.12;

	} else

	if(soiltype == 4){
	/* silt loam */
	tmp.silt = 0.65;
	tmp.clay = 0.15;
	tmp.sand = 0.20;
	tmp.air_entry = -2.1;
	tmp.b = 4.7;
	tmp.Ks = 1.9e-4;
	tmp.satur = 0.59;
	tmp.fieldc = 0.33;
	tmp.wiltp = 0.13;

	} else

	if(soiltype == 5){
	/* sandy clay loam */
	tmp.silt = 0.13;
	tmp.clay = 0.27;
	tmp.sand = 0.60;
	tmp.air_entry = -2.8;
	tmp.b = 4;
	tmp.Ks = 1.2e-4;
	tmp.satur = 0.48;
	tmp.fieldc = 0.26;
	tmp.wiltp = 0.15;

	} else

	if(soiltype == 6){
	/* clay loam */
	tmp.silt = 0.34;
	tmp.clay = 0.34;
        tmp.sand = 0.32;
	tmp.air_entry = -2.6;
	tmp.b = 5.2;
	tmp.Ks = 6.4e-5;
	tmp.satur = 0.52;
	tmp.fieldc = 0.32;
	tmp.wiltp = 0.20;

	} else

	if(soiltype == 7){
	/* silty clay loam */
	tmp.silt = 0.58;
	tmp.clay = 0.33;
	tmp.sand = 0.09;
	tmp.air_entry = -3.3;
	tmp.b = 6.6;
	tmp.Ks = 4.2e-5;
	tmp.satur = 0.52;
	tmp.fieldc = 0.37;
	tmp.wiltp = 0.21; /* Correction from the book from here http://www.public.iastate.edu/~bkh/teaching/505/norman_book_corrections.pdf */

	} else

	if(soiltype == 8){
	/* sandy clay */
	tmp.silt = 0.07;
	tmp.clay = 0.40;
	tmp.sand = 0.53;
	tmp.air_entry = -2.9;
	tmp.b = 6;
	tmp.Ks = 3.3e-5;
	tmp.satur = 0.51;
	tmp.fieldc = 0.34;
	tmp.wiltp = 0.24;

	} else

	if(soiltype == 9){
	/* silty clay */
	tmp.silt = 0.45;
	tmp.clay = 0.45;
	tmp.sand = 0.10;
	tmp.air_entry = -3.4;
	tmp.b = 7.9;
	tmp.Ks = 2.5e-5;
	tmp.satur = 0.52;
	tmp.fieldc = 0.39;
	tmp.wiltp = 0.25;

	} else

	if(soiltype == 10){
	/* clay */
	tmp.silt = 0.20;
	tmp.clay = 0.60;
	tmp.sand = 0.20;
	tmp.air_entry = -3.7;
	tmp.b = 7.6;
	tmp.Ks = 1.7e-5;
	tmp.satur = 0.53;
	tmp.fieldc = 0.4;
	tmp.wiltp = 0.27;

	}

	return(tmp);

}
