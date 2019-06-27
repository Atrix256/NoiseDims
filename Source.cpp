
// Settings


#define DETERMINISTIC() true  // if true, will use the seed below for everything, else will randomly generate a seed.

#define DETERMINISTIC_SEED() unsigned(783104853), unsigned(4213684301), unsigned(3526061164), unsigned(614346169), unsigned(478811579), unsigned(2044310268), unsigned(3671768129), unsigned(206439072)



////////////////////////////////////////////////////////////////

#define _CRT_SECURE_NO_WARNINGS
#include <random>
#include <stdint.h>

#define __USE_SQUARE_BRACKETS_FOR_ELEMENT_ACCESS_OPERATOR
#include "simple_fft/fft_settings.h"
#include "simple_fft/fft.h"

inline std::seed_seq& GetRNGSeed()
{
#if DETERMINISTIC()
    static std::seed_seq fullSeed{ DETERMINISTIC_SEED() };
#else
    static std::random_device rd;
    static std::seed_seq fullSeed{ rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd() };
#endif
    return fullSeed;
}

int RollDice(std::mt19937& rng, int sides)
{
    std::uniform_int_distribution<int> dist(0, sides-1);
    return dist(rng);
}

void DFT(const std::vector<float>& imageSrc, std::vector<float>& imageDest)
{
    // convert the source image to float and store it as complex so it can be DFTd
    size_t width = imageSrc.size();
    std::vector<complex_type> complexImageIn(width);
    for (size_t index = 0; index < width; ++index)
        complexImageIn[index] = imageSrc[index];

    // DFT the image to get frequency of the samples
    const char* error = nullptr;
    std::vector<complex_type> complexImageOut(width);
    simple_fft::FFT(complexImageIn, complexImageOut, width, error);

    // TODO: for some reason, the DC is huge. i'm not sure why...
    complexImageOut[0] = 0.0f;

    // get the magnitudes and max magnitude
    std::vector<float> magnitudes;
    float maxMag = 0.0f;
    {
        magnitudes.resize(width, 0.0f);
        float* dest = magnitudes.data();
        for (size_t x = 0; x < width; ++x)
        {
            size_t srcX = (x + width / 2) % width;

            const complex_type& c = complexImageOut[srcX];
            float mag = float(sqrt(c.real()*c.real() + c.imag()*c.imag()));
            maxMag = std::max(mag, maxMag);
            *dest = mag;
            ++dest;
        }
    }

    // normalize the magnitudes and convert it back to a type T image
    //const float c = 1.0f / log(1.0f / 255.0f + maxMag);
    {
        imageDest.resize(width);
        const float* src = magnitudes.data();
        float* dest = imageDest.data();
        for (size_t x = 0; x < width; ++x)
        {
            //float normalized = c * log(1.0f / 255.0f + *src);
            float normalized = *src / maxMag;
            *dest = normalized;

            ++src;
            ++dest;
        }
    }
}

void TestData(const char* baseFileName, const std::vector<float>& rolls, float mean, float sigma)
{
    // get min / max values
    float minValue = rolls[0];
    float maxValue = rolls[0];
    for (float f : rolls)
    {
        minValue = std::min(minValue, f);
        maxValue = std::max(maxValue, f);
    }

    // make histogram
    int numRolls = int(rolls.size());
    std::vector<int> histogram(101, 0);
    for (int index = 0; index < numRolls; ++index)
    {
        size_t percent = size_t(100.0f * (rolls[index] - minValue) / (maxValue - minValue));
        histogram[percent]++;
    }

    // write histogram
    {
        char fileName[256];
        sprintf(fileName, "%s_%i_%i_%i_histogram.csv", baseFileName, int(mean), int(sigma), numRolls);
        FILE* file = nullptr;
        fopen_s(&file, fileName, "w+t");
        for (int index = 0; index < histogram.size(); ++index)
            fprintf(file, "\"%i\",\"%i\"\n", index, histogram[index]);
        fclose(file);
    }

    // write DFT
    {
        std::vector<float> rollsFloatMag;
        DFT(rolls, rollsFloatMag);

        char fileName[256];
        sprintf(fileName, "%s_%i_%i_%i_DFTMag.csv", baseFileName, int(mean), int(sigma), numRolls);
        FILE* file = nullptr;
        fopen_s(&file, fileName, "w+t");
        for (int index = 0; index < int(numRolls); ++index)
            fprintf(file, "\"%i\",\"%f\"\n", int(index) - int(numRolls) / 2, rollsFloatMag[index]);
        fclose(file);
    }

    // write actual values
    {
        char fileName[256];
        sprintf(fileName, "%s_%i_%i_%i.txt", baseFileName, int(mean), int(sigma), numRolls);
        FILE* file = nullptr;
        fopen_s(&file, fileName, "w+t");
        for (float f : rolls)
            fprintf(file, "%f\n", f);
        fclose(file);
    }
}

void TestData(const char* baseFileName, const std::vector<int>& rolls, int sides, int minValue, int maxValue)
{
    // make float data and histogram
    int numRolls = int(rolls.size());
    std::vector<float> rollsFloat(numRolls);
    std::vector<int> histogram(maxValue - minValue + 1, 0);
    for (int index = 0; index < numRolls; ++index)
    {
        rollsFloat[index] = float(rolls[index]);
        histogram[rolls[index] - minValue]++;
    }

    // write histogram
    {
        char fileName[256];
        sprintf(fileName, "%s_%i_%i_histogram.csv", baseFileName, sides, numRolls);
        FILE* file = nullptr;
        fopen_s(&file, fileName, "w+t");
        for (int index = 0; index < histogram.size(); ++index)
            fprintf(file, "\"%i\",\"%i\"\n", index + minValue, histogram[index]);
        fclose(file);
    }

    // write DFT
    {
        std::vector<float> rollsFloatMag;
        DFT(rollsFloat, rollsFloatMag);

        char fileName[256];
        sprintf(fileName, "%s_%i_%i_DFTMag.csv", baseFileName, sides, numRolls);
        FILE* file = nullptr;
        fopen_s(&file, fileName, "w+t");
        for (int index = 0; index < int(numRolls); ++index)
            fprintf(file, "\"%i\",\"%f\"\n", int(index) - int(numRolls) / 2, rollsFloatMag[index]);
        fclose(file);
    }

    // write actual values
    {
        char fileName[256];
        sprintf(fileName, "%s_%i_%i.txt", baseFileName, sides, numRolls);
        FILE* file = nullptr;
        fopen_s(&file, fileName, "w+t");
        for (int i : rolls)
            fprintf(file, "%i\n", i);
        fclose(file);
    }
}

void RollDiceUncorrelated(int sides, int numRolls)
{
    std::mt19937 rng(GetRNGSeed());

    // do dice rolls
    std::vector<int> rolls(numRolls, 0);
    for (int index = 0; index < numRolls; ++index)
        rolls[index] = RollDice(rng, sides);

    // do tests
    TestData("out/uncorrelated", rolls, sides, 0, 5);
}

void GaussianUncorrelated(float mean, float sigma, int numRolls)
{
    std::mt19937 rng(GetRNGSeed());
    std::normal_distribution<float> dist(mean, sigma);

    // do dice rolls
    std::vector<float> rolls(numRolls, 0);
    for (int index = 0; index < numRolls; ++index)
        rolls[index] = dist(rng);

    // do tests
    TestData("out/gaussianuncorrelated", rolls, mean, sigma);
}


void RollDiceAdditiveUncorrelated(int sides, int numDice, int numRolls)
{
    std::mt19937 rng(GetRNGSeed());

    // do dice rolls
    std::vector<int> rolls(numRolls, 0);
    for (int index = 0; index < numRolls; ++index)
    {
        for (int index2 = 0; index2 < numDice; ++index2)
            rolls[index] += RollDice(rng, sides);
    }

    // do tests
    char fileName[256];
    sprintf(fileName, "out/adduncorrelated%i", numDice);
    TestData(fileName, rolls, sides, 0, numDice * 5);
}

void RollDiceAdditiveCorrelated(int sides, int numDice, int numRolls)
{
    std::mt19937 rng(GetRNGSeed());

    // prime the dice
    std::vector<int> currentDice(numDice);
    for (int& d : currentDice)
        d = RollDice(rng, sides);

    // do dice rolls
    std::vector<int> rolls(numRolls, 0);
    for (int index = 0; index < numRolls; ++index)
    {
        currentDice[index % numDice] = RollDice(rng, sides);

        for (int d : currentDice)
            rolls[index] += d;
    }

    // do tests
    char fileName[256];
    sprintf(fileName, "out/addcorrelated%i", numDice);
    TestData(fileName, rolls, sides, 0, numDice * 5);
}

void RollDiceSubtractUncorrelated(int sides, int numRolls)
{
    int numDice = 2;

    std::mt19937 rng(GetRNGSeed());

    // do dice rolls
    std::vector<int> currentDice(numDice);
    std::vector<int> rolls(numRolls, 0);
    for (int index = 0; index < numRolls; ++index)
    {
        for (int index2 = 0; index2 < numDice; ++index2)
            currentDice[index2] = RollDice(rng, sides);

        rolls[index] = currentDice[index % numDice] - currentDice[(index + 1) % numDice];
    }

    // do tests
    TestData("out/subuncorrelated2", rolls, sides, -(sides - 1), sides - 1);
}

void RollDiceSubtractCorrelated(int sides, int numRolls)
{
    int numDice = 2;

    std::mt19937 rng(GetRNGSeed());

    // prime the dice
    std::vector<int> currentDice(numDice);
    for (int& d : currentDice)
        d = RollDice(rng, sides);

    // do dice rolls
    std::vector<int> rolls(numRolls, 0);
    for (int index = 0; index < numRolls; ++index)
    {
        currentDice[index % numDice] = RollDice(rng, sides);

        rolls[index] = currentDice[index % numDice] - currentDice[(index+1) % numDice];
    }

    // do tests
    TestData("out/subcorrelated2", rolls, sides, -(sides-1), sides-1);
}

int main(int argc, char** argv)
{
    // white noise drawn from rectangular distribution
    printf("RollDiceUncorrelated\n\n");
    RollDiceUncorrelated(6, 16);
    RollDiceUncorrelated(6, 128);
    RollDiceUncorrelated(6, 1024);
    RollDiceUncorrelated(6, 1024*1024);

    // white noise drawn from a gaussian distribution
    printf("GaussianUncorrelated\n\n");
    GaussianUncorrelated(0.0f, 1.0f, 16);
    GaussianUncorrelated(0.0f, 1.0f, 128);
    GaussianUncorrelated(0.0f, 1.0f, 1024);
    GaussianUncorrelated(0.0f, 1.0f, 1024 * 1024);

    // more white noise drawn from a gaussian distribution
    printf("GaussianUncorrelated\n\n");
    GaussianUncorrelated(1.0f, 3.0f, 16);
    GaussianUncorrelated(1.0f, 3.0f, 128);
    GaussianUncorrelated(1.0f, 3.0f, 1024);
    GaussianUncorrelated(1.0f, 3.0f, 1024 * 1024);

    // white noise drawn from a triangular distribution
    printf("RollDiceAdditiveUncorrelated\n\n");
    RollDiceAdditiveUncorrelated(6, 2, 16);
    RollDiceAdditiveUncorrelated(6, 2, 128);
    RollDiceAdditiveUncorrelated(6, 2, 1024);
    RollDiceAdditiveUncorrelated(6, 2, 1024 * 1024);

    // more white noise drawn from a triangular distribution
    printf("RollDiceAdditiveUncorrelated\n\n");
    RollDiceAdditiveUncorrelated(6, 3, 16);
    RollDiceAdditiveUncorrelated(6, 3, 128);
    RollDiceAdditiveUncorrelated(6, 3, 1024);
    RollDiceAdditiveUncorrelated(6, 3, 1024 * 1024);

    // red noise drawn from a triangular distribution
    printf("RollDiceAdditiveCorrelated\n\n");
    RollDiceAdditiveCorrelated(6, 2, 16);
    RollDiceAdditiveCorrelated(6, 2, 128);
    RollDiceAdditiveCorrelated(6, 2, 1024);
    RollDiceAdditiveCorrelated(6, 2, 1024*1024);

    // more red noise drawn from a triangular distribution
    printf("RollDiceAdditiveCorrelated\n\n");
    RollDiceAdditiveCorrelated(6, 3, 16);
    RollDiceAdditiveCorrelated(6, 3, 128);
    RollDiceAdditiveCorrelated(6, 3, 1024);
    RollDiceAdditiveCorrelated(6, 3, 1024 * 1024);

    // blue noise drawn from a triangular distribution
    printf("RollDiceSubtractCorrelated\n\n");
    RollDiceSubtractCorrelated(6, 16);
    RollDiceSubtractCorrelated(6, 128);
    RollDiceSubtractCorrelated(6, 1024);
    RollDiceSubtractCorrelated(6, 1024 * 1024);

    // white noise drawn from a triangular distribution
    printf("RollDiceSubtractUncorrelated\n\n");
    RollDiceSubtractUncorrelated(6, 16);
    RollDiceSubtractUncorrelated(6, 128);
    RollDiceSubtractUncorrelated(6, 1024);
    RollDiceSubtractUncorrelated(6, 1024 * 1024);

    return 0;
}


/*

TODO:
- the experiments from https://www.digido.com/ufaqs/dither-noise-probability-density-explained/
- particularly that drawing from gaussian can be white noise?!


? how would you have more than 2 dice in the RollDiceSubtract?

*/