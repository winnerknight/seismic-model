
/// Stack seismic Forward Modeling With Wave Equation Phase Shift //
// Teaching Programme //
// Writer: XiongGaojun //
// Company: Chengdu University of Technology //
// Date : 2012-05-12 //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Part 1 Preprocedure //
// Nane: PsFrwrdMdlParameter.h //
// Position: Located Inside the same Folder of orther Programmes //
// Function: //
// (1) Various Function Store Include //
// (2) Functions Request: Individual function to be Called must be Request here //
// (3) Constant Symbol Define //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//PsFrwrdMdlParameter.h
// 1. Include Function Store
// 2. Functinns Request
// 3. Parameters set
//---1.Include---------------
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
//---2.Function Request------
int Po2Judge(int); //Judge Nt or Nx is Power of 2
int Absorb(); // Forming Absorbing Boundary
int Rflct(); //Forming Reflect Structure Model
int Vlcty(); //Forming Velociy Model
int kkfft(float*, float*, int, int); // FFT or IFFT
int WvFld0(); //Wave Field Initialization
int exp_ikzDz(float*, int, float, int, float, float);//Compute PhaseShift Data
int PsFrwd(); //Wave Field PhaseShift
//---3. Parameters set---------
#define Nx 128 // Trace Number
#define Nt 256 // Record Number
#define Nz 100 // Depth Number
#define Labs 15 // Length Of Boundary Absorbing
#define Dx 20. // Trace Interval
#define Dt 0.004 // Record Interval
#define Dz 20 // Depth Interval
#define Pai 3.14
//------------------------------
//#include"PsFrwrdMdlParameter.h"
/////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
// Part 2 Main Programme //
// Name: void main( ) //
// Function: Call Functions In the sequence of formula //
///////////////////////////////////////////////////////////////////////////////////////////////
void main()
{
	if (Po2Judge(Nt) != 1) { printf("Nt=%d is not the Power of 2\n", Nt); exit(0); }
	if (Po2Judge(Nx) != 1) { printf("Nx=%d is not the Power of 2\n", Nx); exit(0); }
	if (Absorb() != 1) { printf("Absorb is error\n"); exit(0); }
	if (Rflct() != 1) { printf("Rflction is error\n"); exit(0); }
	if (Vlcty() != 1) { printf("Vlcty is error\n"); exit(0); }
	if (WvFld0() != 1) { printf("WvFld is error\n"); exit(0); }
	if (PsFrwd() != 1) { printf("PsFwrd is error\n"); exit(0); }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Part 3 Functions Programme //
// Function: Individual Function programme Called By Main Programme //
// or Call by other Function Programme //

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//function01: Judge the power of 2
int Po2Judge(int N)
{
	int k = 0;
	long Ln = 0;
	for (k = 0; N - Ln > 0; k++)
	{
		Ln = (long)pow(2, k);
	}
	Ln = (long)pow(2, k - 1);
	if (fabs(Ln - N) >= 1)return(0);
	return(1);
}
////////////////////////////////////////////////////////////////////////////
//function02: Absorb Boundery
int Absorb()
{
	FILE* fp_Abs;
	int Ix;
	float Abs[Nx];
	if ((fp_Abs = fopen("Absb.dat", "wb")) == NULL)printf("Connot open file ""Absb""");
	for (Ix = 0; Ix < Nx; Ix++)
	{
		Abs[Ix] = 1;//1.Absorb Boundary Initializing £¿
	}
	for (Ix = 0; Ix < Labs; Ix++)
	{
		Abs[Ix] = sqrt(sin((Pai / 2) * (Ix / (Labs - 1))));//2.Absorb Boundary Compute£¿
		Abs[Nx - Ix - 1] = Abs[Ix];
	}
	for (Ix = 0; Ix < Nx; Ix++)
	{
		fwrite(&Abs[Ix], sizeof(Abs[Ix]), Nx, fp_Abs);//3.Byte Number of individual Data £¿
	}
	fclose(fp_Abs);
	return(1);
}
///////////////////////////////////////////////////////
//function03: Form Reflect Structure Model
int Rflct()
{
	FILE* fp_Rflct;
	int Ix, Iz;
	float Rflct[Nz];
	if ((fp_Rflct = fopen("Rflct.dat", "wb")) == NULL)printf("Connot open file ""Reflection""");
	for (Ix = 0; Ix < Nx; Ix++)
	{
		for (Iz = 0; Iz < Nz; Iz++)
		{
			Rflct[Iz] = 0.;
			if (Iz == Nz / 5 && fmod(Ix + 1, 2) == 0.) Rflct[Iz] = 1.;
			if (Iz == 2 * Nz / 5 && fmod(Ix + 1, 8) == 0.) Rflct[Iz] = 3.;
			if (Iz == 3 * Nz / 5 && fmod(Ix + 1, 16) == 0.) Rflct[Iz] = 5.;
			if (Iz == 4 * Nz / 5 && fmod(Ix + 1, 32) == 0.) Rflct[Iz] = 7.;
			fwrite(&Rflct[Iz], sizeof(Rflct[Iz]), 1, fp_Rflct);//4.Byte Number of individual Data£¿
		}
	}
	fclose(fp_Rflct);
	return(1);
}
////////////////////////////////////////////////////////
//function04: Form Velocity Model
int Vlcty()
{
	FILE * fp_Vlcty;
	int Ix, Iz;
	float Vlcty[Nz];
	if ((fp_Vlcty = fopen("Vlcty.dat", "wb")) == NULL)printf("Connot open file ""Vlcty""");
	for (Ix = 0; Ix < Nx; Ix++)
	{
		for (Iz = 0; Iz < (int)(3 * Nz / 4); Iz++)
		{
			Vlcty[Iz] = 5000.;
			fwrite(&Vlcty[Iz], sizeof(Vlcty[Iz]), 1, fp_Vlcty);//5.Pointer of Individual data File£¿
		}
		for (Iz = (int)(3 * Nz / 4); Iz < Nz; Iz++)
		{
			Vlcty[Iz] = 5500.;
			fwrite(&Vlcty[Iz], sizeof(Vlcty[Iz]), 1, fp_Vlcty);//6.Byte Number of individual Data£¿
		}
	}
	fclose(fp_Vlcty);
	return(1);
}
////////////////////////////////////////////////////////////////
//function05: Form Inition Wave Field: ReflectCoefficiency //
int WvFld0()
{
	// Function05.1: kkfft(Wfld0r, Wfld0i,Nt,0)
	FILE* fp_Rflct, * fp_Wfld0r, * fp_Wfld0i;
	int Ix, Iz, It;
	float Rflct[Nz];
	float Wfld0r[Nt], Wfld0i[Nt];
	//
	if ((fp_Wfld0r = fopen("Wfld0r.dat", "wb")) == NULL)printf("Connot open Wfld0r.dat");
	if ((fp_Wfld0i = fopen("Wfld0i.dat", "wb")) == NULL)printf("Connot open Wfld0i.dat");
	if ((fp_Rflct = fopen("Rflct.dat", "rb")) == NULL)printf("Connot open Rflct.dat");
	for (Ix = 0; Ix < Nx; Ix++)
	{
		printf("Wavefield0_FFT: Ix=%d\n", Ix);
		for (Iz = 0; Iz < Nz; Iz++)
		{
			fread(&Rflct[Iz], sizeof(Rflct[Iz]), 1, fp_Rflct);
			for (It = 0; It < Nt; It++)
			{
				Wfld0r[It] = 0;//7.Initial Wave field Initializing: Real Part ?
				Wfld0i[It] = 0;//8.Initial Wave field Initializing: Imagine Part ?
			}
			Wfld0r[0] = Rflct[Iz];//9.Initial Wave field Initializing£¿
			if (kkfft(Wfld0r, Wfld0i, Nt, 0) != 1)//10.FFT or Inver FFT ?
			{
				printf("FFT is error");
				exit(0);
			}
			for (It = 0; It < Nt / 2 + 1; It++)//
			{
				fwrite(&Wfld0r[It], sizeof(Wfld0r[It]), 1, fp_Wfld0r);//11.Byte Number of individual Data£¿
				fwrite(&Wfld0i[It], sizeof(Wfld0i[It]), 1, fp_Wfld0i);//12.Byte Number of individual Data£¿
			}
		}
	}
	fclose(fp_Rflct);
	fclose(fp_Wfld0r);
	fclose(fp_Wfld0i);
	return(1);

}
///////////////////////////////////////////////////////////////
//function06: PhaseShift Forward Modeling
int PsFrwd()
{
	//function06.1: PhaseShift calculate
	//function06.2: Wave Field Fourier Transform From Frequency To Time Dpmain
	int PhaseShift();// Request Function:PhaseShift
	int Frqcy2Time();// Request Function:Wave Field Transform From Frequency to Time Domain
	//////
	if (PhaseShift() != 1) { printf("PhaseShift is error\n"); exit(0); }// Call Function
	if (Frqcy2Time() != 1) { printf("Frqcy2Time is error\n"); exit(0); }// Call Function
	return(1);
}
////////////////////////////////////////////////////////////////
//function06.1: PhaseShift calculate
int PhaseShift()
{
	//1 Preprocedure
	//1.1 Functions Request
	//function06.1.1: ReadVlctyAbsb(float *,float *);// Read in Velocity and Absorb Data
	//function06.1.2: FrmNewWfld(FILE *,FILE *,float *,float *,float *,int,int);//Forming New Wave Field
	//function06.1.3: MoveOneDz(float *,float *,float,float,float,int);//Wave Field Extrapolat One Detph Step
	// 
	//int ReadVlctyAbsb(float*, float*);
	int ReadVlctyAbsb(float Vlcty[], float Absb[]);
	//int FrmNewWfld(FILE*, FILE*, float*, float*, float*, int, int);
	int FrmNewWfld(FILE * fp_Wfld0r, FILE * fp_Wfld0i, float Wfldr[], float Wfldi[], float Absb[], int Iz, int Iw);
	//int MoveOneDz(float*, float*, float, float, float, int);
	int MoveOneDz(float Wfldr[], float Wfldi[], float Vz, float Dkx, float Dw, int Iw);
	// 1.2 Define Varibles
	FILE* fp_Wfldr, * fp_Wfldi, * fp_Wfld0r, * fp_Wfld0i;
	int Ix, Iz, Iw, Nw = Nt;
	float Vlcty[Nz], Absb[Nx];
	float Wfldr[Nx], Wfldi[Nx];
	float Kxmax, Dkx, Wmax, Dw;
	// 1.3 Compute out Dw and Dkx
	Dkx = (Pai / Dx) / Nx;//13.Sample Invterval of Wave Number?
	Dw = (Pai / Dt) / Nt;//14.Sample Invterval of Frequency?
	// 2.Read in Velocity and Absorbing Boundary
	if (ReadVlctyAbsb(Vlcty, Absb) != 1) { printf("ReadVlctyAbsb Error"); exit(0); }
	// 3.Open Initial Wave Field File and Current Wave Field File using In Wave Fied Extrapolating
	if ((fp_Wfld0r = fopen("Wfld0r.dat", "rb")) == NULL) printf("Connot open Wfld0r.dat");
	if ((fp_Wfld0i = fopen("Wfld0i.dat", "rb")) == NULL) printf("Connot open Wfld0i.dat");
	if ((fp_Wfldr = fopen("Wfldr.dat", "wb")) == NULL) printf("Connot open Wfldr.dat");
	if ((fp_Wfldi = fopen("Wfldi.dat", "wb")) == NULL) printf("Connot open Wfldi.dat");
	// 4.Wave Fied Extrapolate With Individual Frequency
	for (Iw = 0; Iw < Nw / 2 + 1; Iw++)
	{
		// 4.1 Initializing WaveField of Individual Frequency
		printf("Iw=%d\n", Iw);
		for (Ix = 0; Ix < Nx; Ix++)
		{
			Wfldr[Ix] = 0;//15.PhaseShift Wave field Initializing: Real Part ?
			Wfldi[Ix] = 0;//16.PhaseShift Wave field Initializing: Imagine Part ?
		}
		// 4.2 Wave Fied Extrapolate From Z=Zmax to Z=Zmin
		for (Iz = Nz - 1; Iz > 0; Iz--)
		{
			//4.2.1 Forming New Wave field for Extrapolating on Space Domain
			if (FrmNewWfld(fp_Wfld0r, fp_Wfld0i, Wfldr, Wfldi, Absb, Iz, Iw) != 1) {printf("FrmNewWfld Error");exit(0);}
			// 4.2.3 Wave field Extrapolating One Depth Step Error");exit(0);}
			if (MoveOneDz(Wfldr, Wfldi, (float)(Vlcty[Iz] / 2.), Dkx, Dw, Iw) != 1) {printf("WfldMoveOneStep Error");exit(0);}
		}
		// 4.3 Storing Wave field Data On Survey Line:Z=Zmin
		for (Ix = 0; Ix < Nx; Ix++)
		{
			fwrite(&Wfldr[Ix], sizeof(Wfldr[Ix]), 1, fp_Wfldr);//17.Byte Number of individual Data£¿
			fwrite(&Wfldi[Ix], sizeof(Wfldi[Ix]), 1, fp_Wfldi);//18.Byte Number of individual Data£¿
		}
	}
	///////////////////////////
	// 5 Close or Remove Data File
	fclose(fp_Wfld0r); 
	remove("Wfld0r.dat"); 
	fclose(fp_Wfld0i); 
	remove("Wfld0i.dat");
	fclose(fp_Wfldr);
	fclose(fp_Wfldi);
	return(1);
}
///////////////////////////////////////////////////////////////////
// Function06.1.1: Read In Velocity Data and Absorb Boundary Data
int ReadVlctyAbsb(float Vlcty[], float Absb[])
{
	FILE* fp_Vlcty, * fp_Absb;
	int Iz, Ix;
	if ((fp_Vlcty = fopen("Vlcty.dat", "rb")) == NULL) printf("Connot open Rflct.dat");
	for (Iz = 0; Iz < Nz; Iz++)
	{
		fread(&Vlcty[Iz], sizeof(Vlcty[Iz]), 1, fp_Vlcty);//19.Byte number of individual Data£¿
	}
	fclose(fp_Vlcty);
	if ((fp_Absb = fopen("Absb.dat", "rb")) == NULL) printf("Connot open Absb.dat");
	for (Ix = 0; Ix < Nx; Ix++)
	{
		fread(&Absb[Ix], sizeof(Absb[Ix]), 1, fp_Absb);//20.Byte number of individual Data£¿
	}
	fclose(fp_Absb);
	remove("Absb.dat");
	return(1);
}
///////////////////////////////////////////////////////////
// 06.1.2: Form New Wave Field
int FrmNewWfld(FILE * fp_Wfld0r, FILE * fp_Wfld0i, float Wfldr[], float Wfldi[], float Absb[], int Iz, int Iw)
{
	//function06.1.2.1 ReadIxIzIwToIwIzIx(...)
	int ReadIxIzIwToIwIzIx(FILE * fp_Wfld0r, FILE * fp_Wfld0i, float Wfld0r[], float Wfld0i[], int Iz, int Iw);
	//
	int Ix, Nw = Nt;
	float Wfld0r[Nx], Wfld0i[Nx];
	// Take out Initial Wave Field Data With Individual Depth
	if (ReadIxIzIwToIwIzIx(fp_Wfld0r, fp_Wfld0i, Wfld0r, Wfld0i, Iz, Iw) != 1) {printf("exp_ikzDz is error");exit(0);}
	// Form New Wave Field
	for (Ix = 0; Ix < Nx; Ix++)
	{
		// Compute Current New Wave Field
		Wfldr[Ix] = Wfld0r[Ix] + Wfldr[Ix];//21.Current New Wave Field Compute: Real Part ?
		Wfldi[Ix] = Wfld0i[Ix] + Wfldi[Ix];//22.Current New Wave Field Compute: Imagine Part ?
		// Boundary Obsorb of New Wave Fied
		Wfldr[Ix] = Wfldr[Ix] * Absb[Ix];//23.Current New Wave Field Absorbing: Real Part ?
		Wfldi[Ix] = Wfldi[Ix] * Absb[Ix];//24.Current New Wave Field Absorbing: Imagine Part ?
	}
	return(1);
}
///////////////////////////////////////////////////////////
// 06.1.2.1: Read Data From (Ix,Iz,Iw)Oder to (Iw,Iz,Ix)Oder
int ReadIxIzIwToIwIzIx(FILE * fp_Wfld0r, FILE * fp_Wfld0i, float Wfld0r[], float Wfld0i[], int Iz, int Iw)
{
	int Ix, Nw = Nt;
	long AddfrmStrt;
	for (Ix = 0; Ix < Nx; Ix++)
	{
		// Data Numbers From start Position To Current Position
		AddfrmStrt = (Ix * Nz + Iz) * (Nt / 2 + 1) + Iw;//25.Data number From Its File Head?
		// Byte Numbers From start Position To Current Position
		fseek(fp_Wfld0r, sizeof(Wfld0r[Ix]) * AddfrmStrt, 0);//26. Byte number of Individual Data ?
		fread(&Wfld0r[Ix], sizeof(Wfld0r[Ix]), 1, fp_Wfld0r); //27. Byte number of Individual Data ?
		fseek(fp_Wfld0i, sizeof(Wfld0i[Ix]) * AddfrmStrt, 0);//28. Byte number of Individual Data ?
		fread(&Wfld0i[Ix], sizeof(Wfld0i[Ix]), 1, fp_Wfld0i); //29. Byte number of Individual Data ?
	}
	return(1);
}
/////////////////////////////////////////////////////////////////////////
// 06.1.3: Extrapolate One Depth Step
int MoveOneDz(float Wfldr[], float Wfldi[], float Vz, float Dkx, float Dw, int Iw)
{
	// 06.1.3.1: exp_ikzDz(kz,Ix,Vz,Iw,Dw,Dkx)
	// 06.1.3.2: kkfft(Wfldr, Wfldi,Nx,j)
	int Ikx, Nkx = Nx;
    float kz[2], Wfld_r, Wfld_i;
	if (kkfft(Wfldr, Wfldi, Nx, 0) != 1) //30. FFT or Inverse FFT ?
	{
		printf("FFT is error");
		exit(0);
	}
	for (Ikx = 0; Ikx < Nkx / 2 + 1; Ikx++) //31.The Loop scope of Storing Wave Field IN Wave Number Domain ?
	{
		// 4.2.3.1 Computing Phaseshift Function
		if (exp_ikzDz(kz, Ikx, Vz, Iw, Dw, Dkx) != 1) { printf("exp_ikzDz is error"); exit(0); }
		// 4.2.3.2 WaveField multiply Phaseshift Function
		// Compute WaveField Phaseshift
		Wfld_r = kz[0] * Wfldr[Ikx] - kz[1] * Wfldi[Ikx];//32. WaveField Phaseshift Computing: Real Part ?
		Wfld_i = kz[1] * Wfldr[Ikx] + kz[0] * Wfldi[Ikx];//33. WaveField Phaseshift Computing: Imagine Part ?
		Wfldr[Ikx] = Wfld_r;
	    Wfldi[Ikx] = Wfld_i;
		if (Ikx != 0 && Ikx != Nkx / 2)//34.Condition of WaveField conjugate?
		{
			Wfld_r = kz[0] * Wfldr[Nkx - Ikx] - kz[1] * Wfldi[Nkx - Ikx];//35. WaveField conjugate: Real Part ?
			Wfld_i = kz[1] * Wfldr[Nkx - Ikx] + kz[0] * Wfldi[Nkx - Ikx];//36. WaveField conjugate: Imagine Part ?
			Wfldr[Nkx - Ikx] = Wfld_r;
		    Wfldi[Nkx - Ikx] = Wfld_i;
		}
	}
	// 4.2.4 IFFT of New Wave field From WaveNumber to Space
	if (kkfft(Wfldr, Wfldi, Nkx, 1) != 1)//37. FFT or Inverse FFT ?
	{
		printf("FFT is error");
		exit(0);
	}
	return(1);
}
////////////////////////////////////////////////////////////////
// 06.1.3.1: Compute out PhaseShift Data
int exp_ikzDz(float eikzdz[], int Ix, float Vc, int Iw, float Dw, float Dkx)
{
	float kz;
	eikzdz[0] = 0.;
	eikzdz[1] = 0.;
	kz = sqrt(((Iw * Dw) / Vc) * ((Iw * Dw) / Vc) - (Ix * Dkx) * (Ix * Dkx));//38. Kz computing?
	if (kz > 0)//39. Condition for Kz?
	{
		eikzdz[0] = (float)cos(kz * Dz);//40. Real Part of Exp(-iKzDz)?
		eikzdz[1] = (float)-sin(kz * Dz);//41. Imagine Part of Exp(-iKzDz)?
	}
	return(1);
}
/////////////////////////////////////////////////////////////////////////////////
// Function03: Wave Field Fourier Transform From Frequency Domain to Time Domain
int Frqcy2Time()
{
	FILE* fp_Wfldr, * fp_Wfldi;
	FILE* fp_Record;
	int Ix, It, Iw, Nw = Nt;
	float Wfldtr[Nt], Wfldti[Nt];
	long AddFrmStrt;
	if ((fp_Wfldr = fopen("Wfldr.dat", "rb")) == NULL) {printf("Connot open Wfldr.dat" );exit(0);}
	if ((fp_Wfldi = fopen("Wfldi.dat", "rb")) == NULL) {printf("Connot open Wfldi.dat" );exit(0);}
	if ((fp_Record = fopen("Record.dat", "wb")) == NULL) {printf("Connot open Record.dat");exit(0);}
	for (Ix = 0; Ix < Nx; Ix++)
	{
		for (Iw = 0; Iw < Nw / 2 + 1; Iw++)//42.The Loop scope of Storing Wave Field IN Frequency Domain ?
		{
			AddFrmStrt = Iw * Nx + Ix;//43.Data Number From File Head ?
			fseek(fp_Wfldr, sizeof(Wfldtr[Iw])* AddFrmStrt, 0);//44.Byte number of Individual Data ?
			fread(&Wfldtr[Iw], sizeof(Wfldtr[Iw]), 1, fp_Wfldr); //45.Byte number of Individual Data ?
			fseek(fp_Wfldi, sizeof(Wfldti[Iw]) * AddFrmStrt, 0);//46.Byte number of Individual Data ?
			fread(&Wfldti[Iw], sizeof(Wfldti[Iw]), 1, fp_Wfldi); //47.Byte number of Individual Data ?
			if (Iw != 0 && Iw != Nw / 2) //48.Condition of WaveField conjugate?
			{
				Wfldtr[Nw - Iw] = Wfldtr[Iw]; //49. WaveField conjugat: Real Part ?
				Wfldti[Nw - Iw] = -Wfldti[Iw]; //50. WaveField conjugat: Imagine Part ?
			}
		}
		if (kkfft(Wfldtr, Wfldti, Nw, 1) != 1) //51. FFT or Inverse FFT?
		{
			printf("FFT is error");
			exit(0);
		}
		for (It = 0; It < Nt; It++)
		{
			fwrite(&Wfldtr[It], sizeof(Wfldtr[It]), 1, fp_Record);//52. Byte Number of Individual Data ?
		}
	}
	fclose(fp_Wfldr); remove("Wfldr.dat");
	fclose(fp_Wfldi); remove("Wfldi.dat");
	fclose(fp_Record);
	return(1);
}
//////////////////////////////////////////////////
// Fourier Transform: FFT(i=0) or IFFT(l=1)
int kkfft(float pr[], float pi[], int n, int l)
{
	int it, m, is, i, j, nv, l0, il = 0;
	float p, q, s, vr, vi, poddr, poddi;
	float fr[4096], fi[4096];
	int k = 0;
	long Ln = 0;
	for (k = 0; n - Ln > 0; k++)
	{
		Ln = (long)pow(2, k);
	}
	k = k - 1;
	for (it = 0; it <= n - 1; it++)
	{
		m = it;
		is = 0;
		for (i = 0; i <= k - 1; i++)
		{
			j = m / 2;
			is = 2 * is + (m - 2 * j);
			m = j;
		}
		fr[it] = pr[is];
		fi[it] = pi[is];
	}
	pr[0] = 1.0;
	pi[0] = 0.0;
	p = 6.283185306 / (1.0 * n);
	pr[1] = (float)cos(p);
	pi[1] = -(float)sin(p);
	if (l != 0)
	pi[1] = -pi[1];
	for (i = 2; i <= n - 1; i++)
	{
		p = pr[i - 1] * pr[1];
		q = pi[i - 1] * pi[1];
		s = (pr[i - 1] + pi[i - 1]) * (pr[1] + pi[1]);
		pr[i] = p - q;
		pi[i] = s - p - q;
	}
	for (it = 0; it <= n - 2; it = it + 2)
	{
		vr = fr[it];
		vi = fi[it];
		fr[it] = vr + fr[it + 1];
		fi[it] = vi + fi[it + 1];
		fr[it + 1] = vr - fr[it + 1];
		fi[it + 1] = vi - fi[it + 1];
	}
	m = n / 2;
	nv = 2;
	for (l0 = k - 2; l0 >= 0; l0--)
	{
		m = m / 2;
		nv = 2 * nv;
		for (it = 0; it <= (m - 1) * nv; it = it + nv)
		{
			for (j = 0; j <= (nv / 2) - 1; j++)
			{
				p = pr[m * j] * fr[it + j + nv / 2];
				q = pi[m * j] * fi[it + j + nv / 2];
				s = pr[m * j] + pi[m * j];
				s = s * (fr[it + j + nv / 2] + fi[it + j + nv / 2]);
				poddr = p - q;
				poddi = s - p - q;
				fr[it + j + nv / 2] = fr[it + j] - poddr;
				fi[it + j + nv / 2] = fi[it + j] - poddi;
				fr[it + j] = fr[it + j] + poddr;
				fi[it + j] = fi[it + j] + poddi;
			}
		}
	}
	if (l != 0)
	{
		for (i = 0; i <= n - 1; i++)
		{
			fr[i] = fr[i] / (1.0 * n);
			fi[i] = fi[i] / (1.0 * n);
		}
	}
	if (il != 0)
	{
		for (i = 0; i <= n - 1; i++)
		{
			pr[i] = sqrt(fr[i] * fr[i] + fi[i] * fi[i]);
			if (fabs(fr[i]) < 0.000001 * fabs(fi[i]))
			{
				if ((fi[i] * fr[i]) > 0) 
					pi[i] = 90.0;
				else
					pi[i] = -90.0;
			}
			else
				pi[i] = atan(fi[i] / fr[i]) * 360.0 / 6.283185306;
		}
	}
	for (i = 0; i < n; i++) 
	{ 
		pr[i] = fr[i]; 
		pi[i] = fi[i];
	}
	return (1);
}