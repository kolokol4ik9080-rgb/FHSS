#include "FreqGrid.h"
#include <vector>
#include <complex>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <boost/lexical_cast.hpp>

#define INP_SAMPL_GUARD 32

CFreqGrid::CFreqGrid() : counter(0)
{
	m_dRelFreq = 0.;
	m_fRelFreq = 0.f;

	m_dCeilFreq = 0;
	m_fCeilFreq = 0.f;

	m_pfc64ToneBuf = nullptr;
	m_pfc64ToneBuf2 = nullptr;
	m_pfc32ToneBuf = nullptr;
	m_pfc32ToneBuf2 = nullptr;
	m_psc16ToneBuf = nullptr;
	m_psc16ToneBuf2 = nullptr;

	m_pfc64MultBuf = nullptr;
	m_pfc32MultBuf = nullptr;
	m_psc16MultBuf = nullptr;

	m_pfc64CarriertBuf = nullptr;
	m_pfc32CarriertBuf = nullptr;
	m_psc16CarriertBuf = nullptr;

	m_iInpSampl = 0;
	m_iInpSamp2 = 0;
	m_iInpSamp3 = 0;

	m_iNumSpec = 0;
	m_iWorkType = 0;

	strVal;
	FreqVec;
	
}

CFreqGrid::~CFreqGrid()
{
	Free();
}

bool CFreqGrid::Init(double dSamplFreq, double dRelFreq, double dCeilFreq, int iNumSpec, int iWorkType, std::vector<double> issValues)
{//ИНИЦИАЛИЗАЦИЯ ПАРАМЕТРОВ
	if (dRelFreq < 0) m_dRelFreq = 1. + dRelFreq;
	else m_dRelFreq = dRelFreq;
	m_fRelFreq = (float)m_dRelFreq;

	if (dCeilFreq < 0) m_dCeilFreq = 1. + dCeilFreq;
	else m_dCeilFreq = dCeilFreq;

	m_fCeilFreq = (float)m_dCeilFreq;

	m_iNumSpec = iNumSpec;

	m_iWorkType = iWorkType;

	FreqVec = issValues;

	for (int i = 0; FreqVec.size() > i; i++) {
		FreqVec[i] = FreqVec[i] / dSamplFreq;
		if (FreqVec[i] < 0) FreqVec[i] = FreqVec[i] + 1; // если число в векторе меньше нуля, к нему прибавляется единица для отображения в спектре
	}
	Free();
	return 1;
}

int CFreqGrid::ApplyFreq(double dRelFreq)
{//ПРИМЕНЕНИЕ ПАРАМЕТРОВ
	if (dRelFreq < 0) m_dRelFreq = 1. + dRelFreq;
	else m_dRelFreq = dRelFreq;
	m_fRelFreq = (float)m_dRelFreq;
	return 1;
}

int CFreqGrid::ApplyCFreq(double dCeilFreq)
{//ПРИМЕНЕНИЕ ПАРАМЕТРОВ
	if (dCeilFreq < 0) m_dCeilFreq = 1. + dCeilFreq;
	else m_dCeilFreq =  dCeilFreq;
	m_fCeilFreq =  (float)m_dCeilFreq;
	return 1;
}

//=========================================================================
int CFreqGrid::WorkShort(unsigned char* pucInp, unsigned char* pucOut, int iInpByteLen)
{//ОБРАБОТКА ТИПА SHORT

	int iInpSamplLen = iInpByteLen / 4;
	int iOutByteLen = iInpByteLen;
	float fTmp = 0.f;
	Ipp16sc* I16pucInp = (Ipp16sc*)pucInp;
	Ipp16sc* I16pucOut = (Ipp16sc*)pucOut;
	for (int k = 0; k < iInpSamplLen; k++) {
		I16pucOut[k].re = 0;
		I16pucOut[k].im = 0;
	}
	double m_dPhase = 0.;
	float m_fPhase = (float)m_dPhase;
	int NumSpec = m_iNumSpec;
	if (m_iWorkType == 1) {// режим работы 1 с размером выходного буфера меньше входного в N раз
		int iInpSamplLen_Carrier = iInpSamplLen / NumSpec;
		if (m_iInpSampl < iInpSamplLen_Carrier) {
			Free();
			m_iInpSampl = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для генерации тонального сигнала
			m_psc16ToneBuf = ippsMalloc_16sc(m_iInpSampl);
			m_iInpSamp2 = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для суммирования поднесущих
			m_psc16MultBuf = ippsMalloc_16sc(m_iInpSamp2);
		}
		for (int i = 0; NumSpec > i; i++) {
			ippsTone_16sc(m_psc16ToneBuf, iInpSamplLen_Carrier, 16384, m_fRelFreq * i, &m_fPhase, ippAlgHintAccurate);// генерация тонального сигнала с заданной частотой
			ippsMul_16sc_Sfs(m_psc16ToneBuf, I16pucInp + iInpSamplLen_Carrier * i, m_psc16MultBuf, iInpSamplLen_Carrier, 14);// домножение входного буфера на несущую частоту
			ippsAdd_16sc_ISfs(m_psc16MultBuf, I16pucOut, iInpSamplLen_Carrier, 0);// запись получившихся значений в промежуточный буфер
		}
		ippsTone_16sc(m_psc16CarriertBuf, iInpSamplLen_Carrier, 16384, m_dCeilFreq, &m_fPhase, ippAlgHintAccurate);// генерация тонального сигнала, позволяющего центровать относительно количества спектров    
		ippsMul_16sc_ISfs(m_psc16CarriertBuf, (Ipp16sc*)pucOut, iInpSamplLen_Carrier, 14);// домножение и запись сгенерированного тонального сигнала на получившийся выше сигнал        
		return iOutByteLen / NumSpec;// возвращение выходного буфера с размером меньше в N раз
	}
	else if (m_iWorkType == 2){// режим работы 2 с размером выходного буфера, равным входному
			if (m_iInpSampl < iInpSamplLen) {
				Free();
				m_iInpSampl = iInpSamplLen + INP_SAMPL_GUARD;// создание буфера для генерации тонального сигнала
				m_psc16ToneBuf = ippsMalloc_16sc(m_iInpSampl);
				m_iInpSamp2 = iInpSamplLen + INP_SAMPL_GUARD;// создание буфера для размножения спектров
				m_psc16MultBuf = ippsMalloc_16sc(m_iInpSamp2);
			}
			for (int i = 0; NumSpec > i; i++) {
				ippsTone_16sc(m_psc16ToneBuf, iInpSamplLen, 16384, m_fRelFreq * i, &m_fPhase, ippAlgHintAccurate);// генерация тонального сигнала с заданной частотой
				ippsMul_16sc_Sfs(m_psc16ToneBuf, I16pucInp, m_psc16MultBuf, iInpSamplLen, 14);// домножение входного буфера на несущую частоту
				ippsAdd_16sc_ISfs(m_psc16MultBuf, I16pucOut, iInpSamplLen, 0);// запись получившихся значений в промежуточный буфер
			}
			ippsTone_16sc(m_psc16ToneBuf, iInpSamplLen, 16384, m_dCeilFreq, &m_fPhase, ippAlgHintAccurate);// генерация тонального сигнала, позволяющего центровать относительно количества спектров 
			ippsMul_16sc_ISfs(m_psc16ToneBuf, I16pucOut, iInpSamplLen, 14);// домножение и запись сгенерированного тонального сигнала на получившийся выше сигнал  
			iOutByteLen = iInpSamplLen * 4;
			return iOutByteLen;// возвращение выходного буфера с размером, равным входному
	}
	else if (m_iWorkType == 3) {//режим работы 3 ППРЧ, размер выходного буфера равен размеру входного
		if (FreqVec.size() == 0) {
			int iInpSamplLen_Carrier = iInpSamplLen / NumSpec;
			if (m_iInpSampl < iInpSamplLen_Carrier) {
				Free();
				m_iInpSampl = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для генерации тонального сигнала
				m_psc16ToneBuf = ippsMalloc_16sc(m_iInpSampl);
				m_iInpSamp2 = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для суммирования поднесущих
				m_psc16MultBuf = ippsMalloc_16sc(m_iInpSamp2);
				m_iInpSamp3 = iInpSamplLen + INP_SAMPL_GUARD;// создание буфера для записи в него части отсчетов входного как вид поднесущей
				m_psc16ToneBuf2 = ippsMalloc_16sc(m_iInpSamp3);
			}
			for (int i = 0; NumSpec > i; i++) {
				Ipp32f newm_fRelfreq = m_fRelFreq * i;// присвоение переменной значение из вектора частот
				ippsTone_16sc(m_psc16ToneBuf, iInpSamplLen_Carrier, 16384, newm_fRelfreq, &m_fPhase, ippAlgHintAccurate);// генерация тонального сигнала с заданной частотой
				ippsMul_16sc_Sfs(m_psc16ToneBuf, I16pucInp + iInpSamplLen_Carrier * i, m_psc16MultBuf, iInpSamplLen_Carrier, 14);// домножение входного буфера на несущую частоту
				ippsCopy_16sc(m_psc16MultBuf + i * NumSpec, I16pucOut + i * iInpSamplLen_Carrier, iInpSamplLen_Carrier);// запись получившихся значений в выходной буфер
			}
			ippsTone_16sc(m_psc16ToneBuf2, iInpSamplLen, 16384, m_dCeilFreq, &m_fPhase, ippAlgHintAccurate);// генерация тонального сигнала, позволяющего центровать относительно количества спектров 
			ippsMul_16sc_ISfs(m_psc16ToneBuf2, I16pucOut, iInpSamplLen, 14);// домножение и запись сгенерированного тонального сигнала на получившийся выше сигнал 
		}
		else {
			int iInpSamplLen_Carrier = iInpSamplLen / FreqVec.size();// задаем размер промежуточного буфера для каждой поднесущей
			if (m_iInpSampl < iInpSamplLen_Carrier) {
				Free();
				m_iInpSampl = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для генерации тонального сигнала
				m_psc16ToneBuf = ippsMalloc_16sc(m_iInpSampl);
				m_iInpSamp2 = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для суммирования поднесущих
				m_psc16MultBuf = ippsMalloc_16sc(m_iInpSamp2);
			}
			for (int i = 0; FreqVec.size() > i; i++) {
				Ipp32f newm_fRelfreq = (float)FreqVec[i];// присвоение переменной значение из вектора частот
				ippsTone_16sc(m_psc16ToneBuf, iInpSamplLen_Carrier, 16384, newm_fRelfreq, &m_fPhase, ippAlgHintAccurate);// генерация тонального сигнала с заданной частотой
				ippsMul_16sc_Sfs(m_psc16ToneBuf, I16pucInp + iInpSamplLen_Carrier * i, m_psc16MultBuf, iInpSamplLen_Carrier, 14);// домножение входного буфера на несущую частоту
				ippsCopy_16sc(m_psc16MultBuf + i * FreqVec.size(), I16pucOut + i * iInpSamplLen_Carrier, iInpSamplLen_Carrier);// запись получившихся значений в выходной буфер
			}
		}
		iOutByteLen = iInpSamplLen * 4;
		return iOutByteLen;// возвращение выходного буфера с размером, равным входному
	}
	else if (m_iWorkType == 4) {
		if (FreqVec.size() == 0) {
			if (counter == NumSpec) counter = 0;
			int iInpSamplLen_Carrier = iInpSamplLen;// задаем размер промежуточного буфера для каждой поднесущей
			if (m_iInpSampl < iInpSamplLen_Carrier) {
				Free();
				m_iInpSampl = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для генерации тонального сигнала
				m_psc16ToneBuf = ippsMalloc_16sc(m_iInpSampl);
				m_iInpSamp2 = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для суммирования поднесущих
				m_psc16MultBuf = ippsMalloc_16sc(m_iInpSamp2);
			}
			Ipp32f newm_fRelfreq = m_fRelFreq * counter;
			ippsTone_16sc(m_psc16ToneBuf, iInpSamplLen_Carrier, 16384, newm_fRelfreq, &m_fPhase, ippAlgHintAccurate);// генерация тонального сигнала с заданной частотой
			ippsMul_16sc_Sfs(m_psc16ToneBuf, I16pucInp, m_psc16MultBuf, iInpSamplLen_Carrier, 14);// домножение входного буфера на несущую частоту
			ippsCopy_16sc(m_psc16MultBuf, I16pucOut, iInpSamplLen_Carrier);// запись получившихся значений в выходной буфер
			ippsTone_16sc(m_psc16ToneBuf, iInpSamplLen, 16384, m_dCeilFreq, &m_fPhase, ippAlgHintAccurate);// генерация тонального сигнала, позволяющего центровать относительно количества спектров 
			ippsMul_16sc_ISfs(m_psc16ToneBuf, I16pucOut, iInpSamplLen, 14);// домножение и запись сгенерированного тонального сигнала на получившийся выше сигнал
		}
		else {
			if (counter == FreqVec.size()) counter = 0;
			int iInpSamplLen_Carrier = iInpSamplLen;// задаем размер промежуточного буфера для каждой поднесущей
			if (m_iInpSampl < iInpSamplLen_Carrier) {
				Free();
				m_iInpSampl = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для генерации тонального сигнала
				m_psc16ToneBuf = ippsMalloc_16sc(m_iInpSampl);
				m_iInpSamp2 = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для суммирования поднесущих
				m_psc16MultBuf = ippsMalloc_16sc(m_iInpSamp2);
			}
			Ipp32f newm_fRelfreq = FreqVec[counter];
			ippsTone_16sc(m_psc16ToneBuf, iInpSamplLen_Carrier, 16384, newm_fRelfreq, &m_fPhase, ippAlgHintAccurate);// генерация тонального сигнала с заданной частотой
			ippsMul_16sc_Sfs(m_psc16ToneBuf, I16pucInp, m_psc16MultBuf, iInpSamplLen_Carrier, 14);// домножение входного буфера на несущую частоту
			ippsCopy_16sc(m_psc16MultBuf, I16pucOut, iInpSamplLen_Carrier);// запись получившихся значений в выходной буфер
		}
		counter++;
		iOutByteLen = iInpSamplLen * 4;
		return iOutByteLen;// возвращение выходного буфера с размером, равным входному
	}

}
// методы WorkFloat и WorkDouble выполняют те же самые операциии, что и в WorkShort для каждого из режимов работы
//=========================================================================
int CFreqGrid::WorkFloat(unsigned char* pucInp, unsigned char* pucOut, int iInpByteLen)
{//ОБРАБОТКА ТИПА FLOAT
	int iInpSamplLen = iInpByteLen / 8;
	int iOutByteLen = iInpByteLen;
	Ipp32fc* I32fcpucInp = (Ipp32fc*)pucInp;
	Ipp32fc* I32fcpucOut = (Ipp32fc*)pucOut;
	for (int k = 0; k< iInpSamplLen; k++){
		I32fcpucOut[k].re = 0;
		I32fcpucOut[k].im = 0;
	}
	double m_dPhase = 0.;
	float m_fPhase = (float)m_dPhase;
	int NumSpec = m_iNumSpec;
	if (m_iWorkType == 1) {
		int iInpSamplLen_Carrier = iInpSamplLen / NumSpec;
		if (m_iInpSampl < iInpSamplLen_Carrier) {
			Free();
			m_iInpSampl = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для генерации тонального сигнала
			m_pfc32ToneBuf = ippsMalloc_32fc(m_iInpSampl);
			m_iInpSamp2 = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для суммирования поднесущих
			m_pfc32MultBuf = ippsMalloc_32fc(m_iInpSamp2);
		}
		for (int i = 0; NumSpec > i; i++) {
			ippsTone_32fc(m_pfc32ToneBuf, iInpSamplLen_Carrier, 1.f, m_fRelFreq * i, &m_fPhase, ippAlgHintAccurate);
			ippsMul_32fc(m_pfc32ToneBuf, I32fcpucInp + iInpSamplLen_Carrier * i, m_pfc32MultBuf, iInpSamplLen_Carrier);// сдвинутый
			ippsAdd_32fc_I(m_pfc32MultBuf, I32fcpucOut, iInpSamplLen_Carrier);
		}
			ippsTone_32fc(m_pfc32ToneBuf, iInpSamplLen_Carrier, 1.f, m_dCeilFreq, &m_fPhase, ippAlgHintAccurate);
			ippsMul_32fc_I(m_pfc32ToneBuf, I32fcpucOut, iInpSamplLen_Carrier);
		return iOutByteLen / NumSpec;
	}
	else if (m_iWorkType == 2) {
		if (m_iInpSampl < iInpSamplLen) {
			Free();
			m_iInpSampl = iInpSamplLen + INP_SAMPL_GUARD;// создание буфера для генерации тонального сигнала
			m_pfc32ToneBuf = ippsMalloc_32fc(m_iInpSampl);
			m_iInpSamp2 = iInpSamplLen + INP_SAMPL_GUARD;// создание буфера для размножения спектров
			m_pfc32MultBuf = ippsMalloc_32fc(m_iInpSamp2);
		}
		for (int i = 0; NumSpec > i; i++){
			ippsTone_32fc(m_pfc32ToneBuf, iInpSamplLen, 1.f, m_fRelFreq * i, &m_fPhase, ippAlgHintAccurate);
			ippsMul_32fc(m_pfc32ToneBuf, I32fcpucInp, m_pfc32MultBuf, iInpSamplLen);// сдвинутый
			ippsAdd_32fc_I(m_pfc32MultBuf, I32fcpucOut, iInpSamplLen);
		}
			ippsTone_32fc(m_pfc32ToneBuf, iInpSamplLen, 1.f, m_dCeilFreq, &m_fPhase, ippAlgHintAccurate);
			ippsMul_32fc_I(m_pfc32ToneBuf, I32fcpucOut, iInpSamplLen);
		iOutByteLen = iInpSamplLen * 8;
		return iOutByteLen;
	}
	else if (m_iWorkType == 3) {
		if (FreqVec.size() == 0) {
			int iInpSamplLen_Carrier = iInpSamplLen / NumSpec;
			if (m_iInpSampl < iInpSamplLen_Carrier) {
				Free();
				m_iInpSampl = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для генерации тонального сигнала
				m_pfc32ToneBuf = ippsMalloc_32fc(m_iInpSampl);
				m_iInpSamp2 = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для суммирования поднесущих
				m_pfc32MultBuf = ippsMalloc_32fc(m_iInpSamp2);
				m_iInpSamp3 = iInpSamplLen + INP_SAMPL_GUARD;// создание буфера для размножения спектров
				m_pfc32ToneBuf2 = ippsMalloc_32fc(m_iInpSamp3);
			}
			for (int i = 0; NumSpec > i; i++) {
				Ipp64f newm_fRelfreq = m_fRelFreq * i;
				ippsTone_32fc(m_pfc32ToneBuf, iInpSamplLen_Carrier, 1.f, newm_fRelfreq, &m_fPhase, ippAlgHintAccurate);
				ippsMul_32fc(m_pfc32ToneBuf, I32fcpucInp + iInpSamplLen_Carrier * i, m_pfc32MultBuf, iInpSamplLen_Carrier);// сдвинутый
				ippsCopy_32fc(m_pfc32MultBuf + i * NumSpec, I32fcpucOut + i * iInpSamplLen_Carrier, iInpSamplLen_Carrier);
			}
			ippsTone_32fc(m_pfc32ToneBuf2, iInpSamplLen, 1.f, m_dCeilFreq, &m_fPhase, ippAlgHintAccurate);
			ippsMul_32fc_I(m_pfc32ToneBuf2, I32fcpucOut, iInpSamplLen);
		}
		else {
			int iInpSamplLen_Carrier = iInpSamplLen / FreqVec.size();
			if (m_iInpSampl < iInpSamplLen_Carrier) {
				Free();
				m_iInpSampl = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для генерации тонального сигнала
				m_pfc32ToneBuf = ippsMalloc_32fc(m_iInpSampl);
				m_iInpSamp2 = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для суммирования поднесущих
				m_pfc32MultBuf = ippsMalloc_32fc(m_iInpSamp2);
			}
			for (int i = 0; FreqVec.size() > i; i++) {
				Ipp64f newm_fRelfreq = (float)FreqVec[i];
				ippsTone_32fc(m_pfc32ToneBuf, iInpSamplLen_Carrier, 1.f, newm_fRelfreq, &m_fPhase, ippAlgHintAccurate);
				ippsMul_32fc(m_pfc32ToneBuf, I32fcpucInp + iInpSamplLen_Carrier * i, m_pfc32MultBuf, iInpSamplLen_Carrier);// сдвинутый
				ippsCopy_32fc(m_pfc32MultBuf + i * FreqVec.size(), I32fcpucOut + i * iInpSamplLen_Carrier, iInpSamplLen_Carrier);
			}
		}
		iOutByteLen = iInpSamplLen * 8;
		return iOutByteLen;
	}
	else if (m_iWorkType == 4) {
		if (FreqVec.size() == 0) {
			if (counter == NumSpec) counter = 0;
			int iInpSamplLen_Carrier = iInpSamplLen;
			if (m_iInpSampl < iInpSamplLen_Carrier) {
				Free();
				m_iInpSampl = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для генерации тонального сигнала
				m_pfc32ToneBuf = ippsMalloc_32fc(m_iInpSampl);
				m_iInpSamp2 = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для суммирования поднесущих
				m_pfc32MultBuf = ippsMalloc_32fc(m_iInpSamp2);
			}
			Ipp64f newm_fRelfreq = m_fRelFreq * counter;
			ippsTone_32fc(m_pfc32ToneBuf, iInpSamplLen_Carrier, 1.f, newm_fRelfreq, &m_fPhase, ippAlgHintAccurate);
			ippsMul_32fc(m_pfc32ToneBuf, I32fcpucInp, m_pfc32MultBuf, iInpSamplLen_Carrier);// сдвинутый
			ippsAdd_32fc_I(m_pfc32MultBuf, I32fcpucOut, iInpSamplLen_Carrier);
			ippsTone_32fc(m_pfc32ToneBuf, iInpSamplLen, 1.f, m_dCeilFreq, &m_fPhase, ippAlgHintAccurate);
			ippsMul_32fc_I(m_pfc32ToneBuf, I32fcpucOut, iInpSamplLen);
		}
		else {
			if (counter == FreqVec.size()) counter = 0;
			int iInpSamplLen_Carrier = iInpSamplLen;
			if (m_iInpSampl < iInpSamplLen_Carrier) {
				Free();
				m_iInpSampl = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для генерации тонального сигнала
				m_pfc32ToneBuf = ippsMalloc_32fc(m_iInpSampl);
				m_iInpSamp2 = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для суммирования поднесущих
				m_pfc32MultBuf = ippsMalloc_32fc(m_iInpSamp2);
			}
			Ipp64f newm_fRelfreq = FreqVec[counter];
			ippsTone_32fc(m_pfc32ToneBuf, iInpSamplLen_Carrier, 1.f, newm_fRelfreq, &m_fPhase, ippAlgHintAccurate);
			ippsMul_32fc(m_pfc32ToneBuf, I32fcpucInp, m_pfc32MultBuf, iInpSamplLen_Carrier);// сдвинутый
			ippsAdd_32fc_I(m_pfc32MultBuf, I32fcpucOut, iInpSamplLen_Carrier);
		}
		counter++;
		iOutByteLen = iInpSamplLen * 8;
		
		return iOutByteLen;
	}
}
//=========================================================================
int CFreqGrid::WorkDouble(unsigned char* pucInp, unsigned char* pucOut, int iInpByteLen)
{//ОБРАБОТКА ТИПА DOUBLE
	Ipp64fc* I64fcpucInp = (Ipp64fc*)pucInp;
	Ipp64fc* I64fcpucOut = (Ipp64fc*)pucOut;
	int iInpSamplLen = iInpByteLen / 16;
	int iOutByteLen = iInpByteLen;
	for (int k = 0; k < iInpSamplLen; k++){
		I64fcpucOut[k].re = 0;
		I64fcpucOut[k].im = 0;
	}
	double m_dPhase = 0.;
	int NumSpec = m_iNumSpec;
	if (m_iWorkType == 1) {

		int iInpSamplLen_Carrier = iInpSamplLen / NumSpec;
		if (m_iInpSampl < iInpSamplLen_Carrier) {
			Free();
			m_iInpSampl = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для генерации тонального сигнала
			m_pfc64ToneBuf = ippsMalloc_64fc(m_iInpSampl);
			m_iInpSamp2 = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для суммирования поднесущих
			m_pfc64MultBuf = ippsMalloc_64fc(m_iInpSamp2);
		}
		for (int i = 0; NumSpec > i; i++) {
			ippsTone_64fc(m_pfc64ToneBuf, iInpSamplLen_Carrier, 1.f, m_fRelFreq * i, &m_dPhase, ippAlgHintAccurate);
			ippsMul_64fc(m_pfc64ToneBuf, I64fcpucInp + iInpSamplLen_Carrier * i, m_pfc64MultBuf, iInpSamplLen_Carrier);// сдвинутый
			ippsAdd_64fc_I(m_pfc64MultBuf, I64fcpucOut, iInpSamplLen_Carrier);
		}
			ippsTone_64fc(m_pfc64ToneBuf, iInpSamplLen_Carrier, 1.f, m_dCeilFreq, &m_dPhase, ippAlgHintAccurate);
			ippsMul_64fc_I(m_pfc64ToneBuf, I64fcpucOut, iInpSamplLen_Carrier);
		return iOutByteLen / NumSpec;
	}
	else if (m_iWorkType == 2){
		if (m_iInpSampl < iInpSamplLen) {
			Free();
			m_iInpSampl = iInpSamplLen + INP_SAMPL_GUARD;// создание буфера для генерации тонального сигнала
			m_pfc64ToneBuf = ippsMalloc_64fc(m_iInpSampl);
			m_iInpSamp2 = iInpSamplLen + INP_SAMPL_GUARD;// создание буфера для размножения спектров
			m_pfc64MultBuf = ippsMalloc_64fc(m_iInpSamp2);
		}
		int NumSpec = m_iNumSpec;
		for (int i = 0; NumSpec > i; i++){
			ippsTone_64fc(m_pfc64ToneBuf, iInpSamplLen, 1.f, m_fRelFreq * i, &m_dPhase, ippAlgHintAccurate);
			ippsMul_64fc(m_pfc64ToneBuf, I64fcpucInp, m_pfc64MultBuf, iInpSamplLen);// сдвинутый
			ippsAdd_64fc_I(m_pfc64MultBuf, (Ipp64fc*)pucOut, iInpSamplLen);
		}
		ippsTone_64fc(m_pfc64ToneBuf, iInpSamplLen, 1.f, m_dCeilFreq, &m_dPhase, ippAlgHintAccurate);
		ippsMul_64fc_I(m_pfc64ToneBuf, (Ipp64fc*)pucOut, iInpSamplLen);
		iOutByteLen = iInpSamplLen * 16;
		return iOutByteLen;
	}
	else if (m_iWorkType == 3) {
		if (FreqVec.size() == 0) {
			int iInpSamplLen_Carrier = iInpSamplLen / NumSpec;
			if (m_iInpSampl < iInpSamplLen_Carrier) {
				Free();
				m_iInpSampl = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для генерации тонального сигнала
				m_pfc64ToneBuf = ippsMalloc_64fc(m_iInpSampl);
				m_iInpSamp2 = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для суммирования поднесущих
				m_pfc64MultBuf = ippsMalloc_64fc(m_iInpSamp2);
				m_iInpSamp3 = iInpSamplLen + INP_SAMPL_GUARD;// создание буфера для суммирования поднесущих
				m_pfc64ToneBuf2 = ippsMalloc_64fc(m_iInpSamp3);
			}
			for (int i = 0; NumSpec > i; i++) {
				Ipp64f newm_fRelfreq = m_fRelFreq * i;
				ippsTone_64fc(m_pfc64ToneBuf, iInpSamplLen_Carrier, 1.f, newm_fRelfreq, &m_dPhase, ippAlgHintAccurate);
				ippsMul_64fc(m_pfc64ToneBuf, I64fcpucInp + iInpSamplLen_Carrier * i, m_pfc64MultBuf, iInpSamplLen_Carrier);// сдвинутый
				ippsCopy_64fc(m_pfc64MultBuf + i * NumSpec, I64fcpucOut + i * iInpSamplLen_Carrier, iInpSamplLen_Carrier);
			}
			ippsTone_64fc(m_pfc64ToneBuf2, iInpSamplLen, 1.f, m_dCeilFreq, &m_dPhase, ippAlgHintAccurate);
			ippsMul_64fc_I(m_pfc64ToneBuf2, (Ipp64fc*)pucOut, iInpSamplLen);
		}
		else {
			int iInpSamplLen_Carrier = iInpSamplLen / FreqVec.size();
			if (m_iInpSampl < iInpSamplLen_Carrier) {
				Free();
				m_iInpSampl = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для генерации тонального сигнала
				m_pfc64ToneBuf = ippsMalloc_64fc(m_iInpSampl);
				m_iInpSamp2 = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для суммирования поднесущих
				m_pfc64MultBuf = ippsMalloc_64fc(m_iInpSamp2);
			}
			for (int i = 0; FreqVec.size() > i; i++) {
				Ipp64f newm_fRelfreq = (float)FreqVec[i];
				ippsTone_64fc(m_pfc64ToneBuf, iInpSamplLen_Carrier, 1.f, newm_fRelfreq, &m_dPhase, ippAlgHintAccurate);
				ippsMul_64fc(m_pfc64ToneBuf, I64fcpucInp + iInpSamplLen_Carrier * i, m_pfc64MultBuf, iInpSamplLen_Carrier);// сдвинутый
				ippsCopy_64fc(m_pfc64MultBuf + i * FreqVec.size(), I64fcpucOut + i * iInpSamplLen_Carrier, iInpSamplLen_Carrier);
			}
		}
		iOutByteLen = iInpSamplLen * 16;
		return iOutByteLen;
	}
	else if (m_iWorkType == 4) {
		if (FreqVec.size() == 0) {
			if (counter == NumSpec) counter = 0;
			int iInpSamplLen_Carrier = iInpSamplLen;
			if (m_iInpSampl < iInpSamplLen_Carrier) {
				Free();
				m_iInpSampl = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для генерации тонального сигнала
				m_pfc64ToneBuf = ippsMalloc_64fc(m_iInpSampl);
				m_iInpSamp2 = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для суммирования поднесущих
				m_pfc64MultBuf = ippsMalloc_64fc(m_iInpSamp2);
			}
			Ipp64f newm_fRelfreq = m_fRelFreq * counter;
			ippsTone_64fc(m_pfc64ToneBuf, iInpSamplLen_Carrier, 1.f, newm_fRelfreq, &m_dPhase, ippAlgHintAccurate);
			ippsMul_64fc(m_pfc64ToneBuf, I64fcpucInp, m_pfc64MultBuf, iInpSamplLen_Carrier);// сдвинутый
			ippsCopy_64fc(m_pfc64MultBuf, I64fcpucOut, iInpSamplLen_Carrier);
			ippsTone_64fc(m_pfc64ToneBuf, iInpSamplLen, 1.f, m_dCeilFreq, &m_dPhase, ippAlgHintAccurate);
			ippsMul_64fc_I(m_pfc64ToneBuf, I64fcpucOut, iInpSamplLen);
		}
		else {
			if (counter == FreqVec.size()) counter = 0;
			int iInpSamplLen_Carrier = iInpSamplLen;
			if (m_iInpSampl < iInpSamplLen_Carrier) {
				Free();
				m_iInpSampl = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для генерации тонального сигнала
				m_pfc64ToneBuf = ippsMalloc_64fc(m_iInpSampl);
				m_iInpSamp2 = iInpSamplLen_Carrier + INP_SAMPL_GUARD;// создание буфера для суммирования поднесущих
				m_pfc64MultBuf = ippsMalloc_64fc(m_iInpSamp2);
			}
			Ipp64f newm_fRelfreq = FreqVec[counter];
			ippsTone_64fc(m_pfc64ToneBuf, iInpSamplLen_Carrier, 1.f, newm_fRelfreq, &m_dPhase, ippAlgHintAccurate);
			ippsMul_64fc(m_pfc64ToneBuf, I64fcpucInp, m_pfc64MultBuf, iInpSamplLen_Carrier);// сдвинутый
			ippsCopy_64fc(m_pfc64MultBuf, I64fcpucOut, iInpSamplLen_Carrier);
		}
		counter++;
		iOutByteLen = iInpSamplLen * 16;
		return iOutByteLen;
	}
}
//=========================================================================
int CFreqGrid::Free()
{
	//ОСВОБОЖДЕНИЕ ПАМЯТИ КЛАССА

	return 1;
}
//=========================================================================