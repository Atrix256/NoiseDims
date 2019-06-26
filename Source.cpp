
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

void TestData(const char* baseFileName, const std::vector<int>& rolls, int sides, int numRolls, int minValue, int maxValue)
{
    // make float data and histogram
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
            fprintf(file, "\"%i\",\"%f\"\n", int(index) - int(numRolls) / 2, rollsFloatMag[index]); // TODO: are these frequencies correct? test & verify!
        // TODO: they definitely aren't correct now that we have a minValue! fix!
        fclose(file);
    }


    // TODO: write out actual values too. Might be useful & fun for some people
}

void RollDiceUncorrelated(int sides, int numRolls)
{
    std::mt19937 rng(GetRNGSeed());

    // do dice rolls
    std::vector<int> rolls(numRolls, 0);
    for (int index = 0; index < numRolls; ++index)
        rolls[index] = RollDice(rng, sides);

    // do tests
    TestData("out/uncorrelated", rolls, sides, numRolls, 0, 5);
}

void RollDiceAdditive(int sides, int numDice, int numRolls)
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
    sprintf(fileName, "out/add%i", numDice);
    TestData(fileName, rolls, sides, numRolls, 0, numDice * 5);
}

void RollDiceSubtract(int sides, int numRolls)
{
    // TODO: how to handle more than 2 dice?
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
    TestData("out/sub2", rolls, sides, numRolls, -5, 5);
}

int main(int argc, char** argv)
{
    // white noise drawn from rectangular distribute
    RollDiceUncorrelated(6, 16);
    RollDiceUncorrelated(6, 128);
    RollDiceUncorrelated(6, 1024);
    RollDiceUncorrelated(6, 1024*1024);

    // TODO: explain. red noise but...
    RollDiceAdditive(6, 2, 16);
    RollDiceAdditive(6, 2, 128);
    RollDiceAdditive(6, 2, 1024);
    RollDiceAdditive(6, 2, 1024*1024);

    RollDiceAdditive(6, 3, 16);
    RollDiceAdditive(6, 3, 128);
    RollDiceAdditive(6, 3, 1024);
    RollDiceAdditive(6, 3, 1024 * 1024);

    // TODO: explain. blue noise but...
    RollDiceSubtract(6, 16);
    RollDiceSubtract(6, 128);
    RollDiceSubtract(6, 1024);
    RollDiceSubtract(6, 1024 * 1024);

    return 0;
}


/*

TODO:
- the experiments from https://www.digido.com/ufaqs/dither-noise-probability-density-explained/
- particularly that drawing from gaussian can be white noise?!
- deal with todos
- is there a way we can zip the csv's to take up less space in repo?

*/