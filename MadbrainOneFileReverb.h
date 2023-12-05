#pragma once

#define MADBRAIN_ONE_FILE_REVERB_FIXED_DELAY_LINES 1

#ifdef MADBRAIN_ONE_FILE_REVERB_FIXED_DELAY_LINES
#include <assert.h>
#else
#include <vector>
#endif

#include <string.h>
#include <math.h>

#define MADBRAIN_ONE_FILE_REVERB_DEFAULT_SAMPLE_RATE 48000

class MadbrainOneFileReverb;

// MadbrainOneFileReverb functions:
// 
// void setSampleRate (float sr) {m_sampleRate  = sr; m_clock = 1;}
// void generate(float *leftInOut, float *rightInOut, int nbSamples)
// void setReverbTime (float v ) // 0 to 1
// void setReverbLevel(float v ) // 0 to 1
// void setStereoWidth(float v ) // 0 to 1
//
// By Hubert "madbrain" Lamontagne in 2022
// Allpass-loop reverb
// Based on Jon Dattorro's clasic paper ( https://ccrma.stanford.edu/~dattorro/ )
//
// I can compose nice music for your game or other media,
// See https://soundcloud.com/madbr for demo songs
// contact me at https://twitter.com/MrMadbrain if interested
//
// MIT license

// longest delay line is 200ms
// 200ms at 48khz needs 9504 samples, rounds up to 16384
#define MADBRAIN_ONE_FILE_REVERB_DEFAULT_DELAY_SIZE 16384  // must always be a power of 2

class MadbrainOneFileReverbDelayLine {
public:
	MadbrainOneFileReverbDelayLine() {
		m_ptr    = 0;
		m_apMem  = 0;
		m_apCoef = 0;
		m_apRamp = 0;
#ifdef MADBRAIN_ONE_FILE_REVERB_FIXED_DELAY_LINES
		memset(&m_data[0], 0, MADBRAIN_ONE_FILE_REVERB_DEFAULT_DELAY_SIZE * sizeof(float));
#else
		m_data.reserve(MADBRAIN_ONE_FILE_REVERB_DEFAULT_DELAY_SIZE);
		m_data.resize(MADBRAIN_ONE_FILE_REVERB_DEFAULT_DELAY_SIZE);
#endif
		m_mask = MADBRAIN_ONE_FILE_REVERB_DEFAULT_DELAY_SIZE - 1;
		m_length = 2;
		m_jump = true;
	}

	inline float output() {
		float dl =  - m_data[(m_ptr - m_length) & m_mask] * m_apCoef + m_data[(m_ptr - m_length - 1) & m_mask];
		m_apMem *= m_apCoef;
		m_apMem += dl;
		return m_apMem;
	}

	inline float outputNoSubsample() {
		float dl =  m_data[(m_ptr - m_length) & m_mask];
		return dl;
	}

	inline void input(float val) {
		m_data[m_ptr & m_mask] = val;
		m_ptr++;

		m_apCoef += m_apRamp;
		while(m_apRamp > 0 && m_apCoef >= 0.2f) {
			m_length++;
			m_apCoef -= 8.f/15;
		}
		while(m_apRamp < 0 && m_apCoef < -1.f / 3.f) {
			m_length--;
			m_apCoef += 8.f/15;
		}
	}

	inline void inputNoSubsample(float val) {
		m_data[m_ptr & m_mask] = val;
		m_ptr++;
	}

	inline void setTimeRamp(float time, float rampTime, bool noSubsample) {
		// check if min clamp should be 1.1 (orig min)
		if(time < 2.001f)
			time = 2.001f;

#ifdef MADBRAIN_ONE_FILE_REVERB_FIXED_DELAY_LINES
		assert(time <= m_mask);
#else
		if(time > m_mask) {
			int oldSize = m_mask + 1;
			while(time > m_mask)
				m_mask += m_mask + 1;
			m_data.reserve(m_mask + 1);
			m_data.resize (m_mask + 1, m_data.back());
			for(int h = 0; h < oldSize; h++)
				m_data[h + m_mask + 1 - oldSize] = m_data[h];
		}
#endif

		int samples  = (int)(time - 1.5f) + 1; // <0 rounds upwards!
		float frac   = time - samples;
		float param  = (frac - 1) / (frac + 1);

		if(m_jump || noSubsample) {
			m_apCoef = param;
			m_length = samples;
			m_apRamp = 0;
			m_jump   = false;
		}
		else {
			float ramp   = param - m_apCoef + (samples - m_length)*(8.f/15);
			ramp /= rampTime;
			m_apRamp = ramp;
		}
	}

	inline void clear() {
#ifdef MADBRAIN_ONE_FILE_REVERB_FIXED_DELAY_LINES
		memset(&m_data[0], 0, MADBRAIN_ONE_FILE_REVERB_DEFAULT_DELAY_SIZE * sizeof(float));
#else
		memset(&m_data[0], 0, m_data.size() * sizeof(float));
#endif
		m_apMem = 0;
	}

	int m_ptr;
	int m_mask;
	int m_length;
	float m_apMem;
	float m_apCoef;
	float m_apRamp;
	bool m_jump;

#ifdef MADBRAIN_ONE_FILE_REVERB_FIXED_DELAY_LINES
	float m_data[MADBRAIN_ONE_FILE_REVERB_DEFAULT_DELAY_SIZE];
#else
	std::vector<float> m_data;
#endif
};


class MadbrainOneFileReverb {
public:

	MadbrainOneFileReverb() {
		static const double lfoInitPhases[] = {
			0.410548654, 0.379876636, 0.6297889578, 0.3081143567, 0.6133103865, 0.7817358636, 0.7612043567, 0.03318636735
		};
		for(int i=0; i<8; i++)
			m_lfoPhase[i] = lfoInitPhases[i];

		m_sampleRate = MADBRAIN_ONE_FILE_REVERB_DEFAULT_SAMPLE_RATE;
		m_clock   = 1;

		m_reverbTime  = 0.65f;
		m_reverbLevel = 0.15f;
		m_stereoWidth = 0.65f;

		m_mix     = 0;
		m_mixR    = 0;
		m_stereo  = 0;
		m_stereoR = 0;

		for(int i=0; i<4; i++) {
			m_filterCoef[i] = 0;
			m_filterMem [i] = 0;
		}
		for(int i=0; i<12; i++) {
			m_delayWet[i] = 0;
		}
		for(int i=0; i<8; i++) {
			m_delayDry[i] = 0;
			m_delayFb [i] = 0;
		}
	}

	inline void setSampleRate (float sr) {m_sampleRate  = sr; m_clock = 1;}
	inline void setReverbTime (float v ) {m_reverbTime  = v; }
	inline void setReverbLevel(float v ) {m_reverbLevel = v; }
	inline void setStereoWidth(float v ) {m_stereoWidth = v; }
	inline void generate(float *leftInOut, float *rightInOut, int nbSamples, int stride = 1) {
		static const double lfoRates[] = {
			0.2233456786543217633, 0.2468753235468365261, 0.0837534234567563417, 0.1566524565321462346,
			0.1937645365732735254, 0.0626735673583567271, 0.1046735673563457619, 0.1362456245756726345
		};
		static const double lfoLevels  [] = {0.9, 0.6, 2.4, 1.2, 0.9, 3.0, 2.1, 1.4 };
		static const double apfTimes   [] = {47.6784, 14.352, 129.6,   52.45, 33.462,  24.567, 113.873, 63.3415 };
		static const double apfCoefs   [] = {0.6, 0.7, 0.4, 0.5, 0.6, 0.7, 0.4, 0.5 };
		static const double delayTimes [] = {21.62, 52.525, 18.42, 87.3454 };
		static const double filterFreqs[] = {8000, 8000, 4000, 100 };
		static const double timeScale = 1.5f;

		for(int i=0; i<nbSamples; i++) {
			m_clock--;
			if(m_clock <= 0) {
				int clockTime = (int)(m_sampleRate / 44100.0 * 32 + 0.5);
				m_clock = clockTime;

				for(int j=0; j<8; j++) {
					float lfo = (float)m_lfoPhase[j];
					lfo = lfo < 0.25f ? lfo*4 : (lfo < 0.75f ? (0.5f-lfo)*4 : (lfo-1)*4);
					m_lfoPhase[j] += lfoRates[j] / m_sampleRate * clockTime;
					if(m_lfoPhase[j] >= 1)
						m_lfoPhase[j] -= 1;
					if(m_lfoPhase[j] >= 1)
						m_lfoPhase[j] = 0;
					double time = m_sampleRate * 0.001f * (apfTimes[j] + lfoLevels[j] * lfo) * timeScale;
					m_delays[j].setTimeRamp((float)time, (float)clockTime, false);
					double absorbCoef = pow(0.001, apfTimes[j] * timeScale * 0.001 / (m_reverbTime + 0.000001));
					m_delayDry[j] = (float)-apfCoefs[j];
					m_delayWet[j] = (float)(1 - apfCoefs[j]*apfCoefs[j]) * (float)absorbCoef;
					m_delayFb [j] = (float)apfCoefs[j] * (float)absorbCoef;
				}
				for(int j=0; j<4; j++) {
					double time = m_sampleRate * 0.001f * delayTimes[j] * timeScale;
					m_delays[j+8].setTimeRamp((float)time, (float)clockTime, true);
					m_delayWet[j+8] = (float)pow(0.001, delayTimes[j] * timeScale * 0.001 / (m_reverbTime + 0.000001));
				}

				float fxMix    = m_reverbLevel;
				float fxStereo = m_stereoWidth * fxMix;
				m_mixR    = (fxMix    - m_mix   ) / clockTime;
				m_stereoR = (fxStereo - m_stereo) / clockTime;

				for(int j=0; j<4; j++) {
					double c = expf(-2 * 3.1415926535f * (float)filterFreqs[j] / m_sampleRate);
					c = c >= 0 && c <= 1 ? c : (c <= 0.5 ? 0.0 : 1.0);
					m_filterCoef[j] = 1 - (float)c;
				}
			}

			float lval = *leftInOut;
			float rval = *rightInOut;
			float inM = ( lval + rval) * 0.5f + 0.000000001f;
			float inS = (-lval + rval) * 0.5f + 0.000000001f;
			float t;
			float o;
			float l = 0;
			float r = 0;

			m_filterMem[2] += (inM            - m_filterMem[2]) * m_filterCoef[2];
			m_filterMem[3] += (m_filterMem[2] - m_filterMem[3]) * m_filterCoef[3];
			float in  = m_filterMem[2] - m_filterMem[3];

			o = m_delays[11].output();
			t = o * m_delayWet[11];

			m_filterMem[0] += (t - m_filterMem[0]) * m_filterCoef[0];
			t = m_filterMem[0];

			l += t;
			t += in;

			o = m_delays[0].output();
			m_delays[0].input(t + o*m_delayFb[0]);
			t = t*m_delayDry[0] + o*m_delayWet[0];

			o = m_delays[1].output();
			m_delays[1].input(t + o*m_delayFb[1]);
			t = t*m_delayDry[1] + o*m_delayWet[1];

			l += t;
			t += in;

			o = m_delays[8].output();
			m_delays[8].input(t);
			t = o * m_delayWet[8];

			r += t;
			t += in;

			o = m_delays[2].output();
			m_delays[2].input(t + o*m_delayFb[2]);
			t = t*m_delayDry[2] + o*m_delayWet[2];

			o = m_delays[3].output();
			m_delays[3].input(t + o*m_delayFb[3]);
			t = t*m_delayDry[3] + o*m_delayWet[3];

			r += t;
			t += in;

			o = m_delays[9].output();
			m_delays[9].input(t);
			t = o * m_delayWet[9];

			m_filterMem[1] += (t - m_filterMem[1]) * m_filterCoef[1];
			t = m_filterMem[1];

			r += t;
			t += in;

			o = m_delays[4].output();
			m_delays[4].input(t + o*m_delayFb[4]);
			t = t*m_delayDry[4] + o*m_delayWet[4];

			o = m_delays[5].output();
			m_delays[5].input(t + o*m_delayFb[5]);
			t = t*m_delayDry[5] + o*m_delayWet[5];

			r += t;
			t += in;

			o = m_delays[10].output();
			m_delays[10].input(t);
			t = o * m_delayWet[10];

			l += t;
			t += in;

			o = m_delays[6].output();
			m_delays[6].input(t + o*m_delayFb[6]);
			t = t*m_delayDry[6] + o*m_delayWet[6];

			o = m_delays[7].output();
			m_delays[7].input(t + o*m_delayFb[7]);
			t = t*m_delayDry[7] + o*m_delayWet[7];

			l += t;
			t += in;

			m_delays[11].input(t);

			float out1 = inM + ( l + r)*m_mix;
			float out2 = inS + (-l + r)*m_stereo;

			*leftInOut  = out1 - out2;
			*rightInOut = out1 + out2;

			leftInOut  += stride;
			rightInOut += stride;

			m_mix    += m_mixR;
			m_stereo += m_stereoR;
		}
	}

	float m_sampleRate;
	int m_clock;

	float m_reverbTime;
	float m_reverbLevel;
	float m_stereoWidth;

	float m_mix;
	float m_mixR;
	float m_stereo;
	float m_stereoR;

	float m_filterCoef[4];
	float m_filterMem [4];

	double m_lfoPhase[8];
	MadbrainOneFileReverbDelayLine m_delays[12];

	float m_delayWet[12];
	float m_delayDry[8];
	float m_delayFb[8];
};


