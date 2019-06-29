
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

void RollDiceAdditiveCorrelated(int sides, int numDice, int numRolls, int numRerolls)
{
    std::mt19937 rng(GetRNGSeed());

    // no funny business boyos
    numRerolls = std::min(numRerolls, numDice);

    // prime the dice
    std::vector<int> currentDice(numDice);
    for (int& d : currentDice)
        d = RollDice(rng, sides);

    int nextToReroll = 0;

    // do dice rolls
    std::vector<int> rolls(numRolls, 0);
    for (int index = 0; index < numRolls; ++index)
    {
        // reroll however many dice we should
        for (int rerollIndex = 0; rerollIndex < numRerolls; ++rerollIndex)
        {
            currentDice[nextToReroll] = RollDice(rng, sides);
            nextToReroll = (nextToReroll + 1) % numDice;
        }

        for (int d : currentDice)
            rolls[index] += d;
    }

    // do tests
    char fileName[256];
    sprintf(fileName, "out/addcorrelated%i_%i", numDice, numRerolls);
    TestData(fileName, rolls, sides, 0, numDice * 5);
}

void GaussianAdditiveCorrelated(int numDice, float mean, float sigma, int numRolls, int numRerolls)
{
    std::mt19937 rng(GetRNGSeed());
    std::normal_distribution<float> dist(mean, sigma);

    // no funny business boyos
    numRerolls = std::min(numRerolls, numDice);

    // prime the dice
    std::vector<float> currentDice(numDice);
    for (float& d : currentDice)
        d = dist(rng);

    int nextToReroll = 0;

    // do dice rolls
    std::vector<float> rolls(numRolls, 0);
    for (int index = 0; index < numRolls; ++index)
    {
        // reroll however many dice we should
        for (int rerollIndex = 0; rerollIndex < numRerolls; ++rerollIndex)
        {
            currentDice[nextToReroll] = dist(rng);
            nextToReroll = (nextToReroll + 1) % numDice;
        }

        for (float d : currentDice)
            rolls[index] += d;
    }

    // do tests
    char fileName[256];
    sprintf(fileName, "out/gaussianaddcorrelated%i_%i", numDice, numRerolls);
    TestData(fileName, rolls, mean, sigma);
}

void RollDiceSubtractUncorrelated(int sides, int numRolls, int numDice)
{
    std::mt19937 rng(GetRNGSeed());

    // do dice rolls
    std::vector<int> currentDice(numDice);
    std::vector<int> rolls(numRolls, 0);
    for (int index = 0; index < numRolls; ++index)
    {
        for (int index2 = 0; index2 < numDice; ++index2)
            currentDice[index2] = RollDice(rng, sides);

        for (int i = 0; i < numDice; ++i)
            rolls[index] += currentDice[i] * ((((index + i) % numDice) % 2) ? -1 : 1);
    }

    // do tests
    char fileName[256];
    sprintf(fileName, "out/subuncorrelated%i", numDice);
    TestData(fileName, rolls, sides, -(sides - 1)*(numDice-1), (sides - 1)*(numDice-1));
}

void RollDiceSubtractCorrelated(int sides, int numRolls, int numDice, int numRerolls)
{
    std::mt19937 rng(GetRNGSeed());

    // no funny business boyos
    numRerolls = std::min(numRerolls, numDice);

    // prime the dice
    std::vector<int> currentDice(numDice);
    for (int& d : currentDice)
        d = RollDice(rng, sides);

    int nextToReroll = 0;

    // do dice rolls
    std::vector<int> rolls(numRolls, 0);
    for (int index = 0; index < numRolls; ++index)
    {
        // reroll however many dice we should
        for (int rerollIndex = 0; rerollIndex < numRerolls; ++rerollIndex)
        {
            currentDice[nextToReroll] = RollDice(rng, sides);
            nextToReroll = (nextToReroll + 1) % numDice;
        }

        for (int i = 0; i < numDice; ++i)
            rolls[index] += currentDice[i] * ((((index + i) % numDice) % 2) ? -1 : 1);
    }

    // do tests
    char fileName[256];
    sprintf(fileName, "out/subcorrelated%i_%i", numDice, numRerolls);
    TestData(fileName, rolls, sides, -(sides - 1)*(numDice - 1), (sides - 1)*(numDice - 1));
}

void GaussianSubtractCorrelated(int numDice, float mean, float sigma, int numRolls, int numRerolls)
{
    std::mt19937 rng(GetRNGSeed());
    std::normal_distribution<float> dist(mean, sigma);

    // no funny business boyos
    numRerolls = std::min(numRerolls, numDice);

    // prime the dice
    std::vector<float> currentDice(numDice);
    for (float& d : currentDice)
        d = dist(rng);

    int nextToReroll = 0;

    // do dice rolls
    std::vector<float> rolls(numRolls, 0);
    for (int index = 0; index < numRolls; ++index)
    {
        // reroll however many dice we should
        for (int rerollIndex = 0; rerollIndex < numRerolls; ++rerollIndex)
        {
            currentDice[nextToReroll] = dist(rng);
            nextToReroll = (nextToReroll + 1) % numDice;
        }

        for (int i = 0; i < numDice; ++i)
            rolls[index] += currentDice[i] * ((((index + i) % numDice) % 2) ? -1 : 1);
    }

    // do tests
    char fileName[512];
    sprintf(fileName, "out/gaussiansubcorrelated%i_%i", numDice, numRerolls);
    TestData(fileName, rolls, mean, sigma);
}

void AddBlueRed(int numRolls)
{
    std::mt19937 rng(GetRNGSeed());

    int numRerolls = 1;
    int sides = 6;
    int numDice = 20;

    std::vector<int> rollsRed(numRolls, 0);
    std::vector<int> rollsBlue(numRolls, 0);

    // red noise
    {
        // prime the dice
        std::vector<int> currentDice(numDice);
        for (int& d : currentDice)
            d = RollDice(rng, sides);

        int nextToReroll = 0;

        // do dice rolls
        for (int index = 0; index < numRolls; ++index)
        {
            // reroll however many dice we should
            for (int rerollIndex = 0; rerollIndex < numRerolls; ++rerollIndex)
            {
                currentDice[nextToReroll] = RollDice(rng, sides);
                nextToReroll = (nextToReroll + 1) % numDice;
            }

            for (int d : currentDice)
                rollsRed[index] += d;
        }
    }

    // blue noise
    {
        // prime the dice
        std::vector<int> currentDice(numDice);
        for (int& d : currentDice)
            d = RollDice(rng, sides);

        int nextToReroll = 0;

        // do dice rolls
        std::vector<int> rolls(numRolls, 0);
        for (int index = 0; index < numRolls; ++index)
        {
            // reroll however many dice we should
            for (int rerollIndex = 0; rerollIndex < numRerolls; ++rerollIndex)
            {
                currentDice[nextToReroll] = RollDice(rng, sides);
                nextToReroll = (nextToReroll + 1) % numDice;
            }

            for (int i = 0; i < numDice; ++i)
                rollsBlue[index] += currentDice[i] * ((((index + i) % numDice) % 2) ? -1 : 1);
        }
    }

    std::vector<int> rolls(numRolls, 0);
    int minValue = 10000;
    int maxValue = -10000;
    for (int index = 0; index < numRolls; ++index)
    {
        rolls[index] = rollsRed[index] + rollsBlue[index];
        minValue = std::min(minValue, rolls[index]);
        maxValue = std::max(maxValue, rolls[index]);
    }

    // do tests
    char fileName[256];
    sprintf(fileName, "out/addbluered%i", numDice);
    TestData(fileName, rolls, 6, minValue, maxValue);
}

void MixBlueRed(int numRolls, float redPercentage)
{
    std::mt19937 rng(GetRNGSeed());

    int numRerolls = 1;
    int sides = 6;
    int numDice = 20;

    std::vector<int> rollsRed(numRolls, 0);
    std::vector<int> rollsBlue(numRolls, 0);

    // red noise
    {
        // prime the dice
        std::vector<int> currentDice(numDice);
        for (int& d : currentDice)
            d = RollDice(rng, sides);

        int nextToReroll = 0;

        // do dice rolls
        for (int index = 0; index < numRolls; ++index)
        {
            // reroll however many dice we should
            for (int rerollIndex = 0; rerollIndex < numRerolls; ++rerollIndex)
            {
                currentDice[nextToReroll] = RollDice(rng, sides);
                nextToReroll = (nextToReroll + 1) % numDice;
            }

            for (int d : currentDice)
                rollsRed[index] += d;
        }
    }

    // blue noise
    {
        // prime the dice
        std::vector<int> currentDice(numDice);
        for (int& d : currentDice)
            d = RollDice(rng, sides);

        int nextToReroll = 0;

        // do dice rolls
        std::vector<int> rolls(numRolls, 0);
        for (int index = 0; index < numRolls; ++index)
        {
            // reroll however many dice we should
            for (int rerollIndex = 0; rerollIndex < numRerolls; ++rerollIndex)
            {
                currentDice[nextToReroll] = RollDice(rng, sides);
                nextToReroll = (nextToReroll + 1) % numDice;
            }

            for (int i = 0; i < numDice; ++i)
                rollsBlue[index] += currentDice[i] * ((((index + i) % numDice) % 2) ? -1 : 1);
        }
    }

    std::vector<int> rolls(numRolls, 0);
    int minValue = 10000;
    int maxValue = -10000;
    for (int index = 0; index < numRolls; ++index)
    {
        rolls[index] = int(float(rollsRed[index]) * 0.5f + float(rollsBlue[index]) * 0.5f);
        minValue = std::min(minValue, rolls[index]);
        maxValue = std::max(maxValue, rolls[index]);
    }

    // do tests
    char fileName[256];
    sprintf(fileName, "out/mixbluered%i_%i", numDice, int(100.0f * redPercentage));
    TestData(fileName, rolls, 6, minValue, maxValue);
}

void MaxBlueRed(int numRolls)
{
    std::mt19937 rng(GetRNGSeed());

    int numRerolls = 1;
    int sides = 6;
    int numDice = 20;

    std::vector<int> rollsRed(numRolls, 0);
    std::vector<int> rollsBlue(numRolls, 0);

    // red noise
    {
        // prime the dice
        std::vector<int> currentDice(numDice);
        for (int& d : currentDice)
            d = RollDice(rng, sides);

        int nextToReroll = 0;

        // do dice rolls
        for (int index = 0; index < numRolls; ++index)
        {
            // reroll however many dice we should
            for (int rerollIndex = 0; rerollIndex < numRerolls; ++rerollIndex)
            {
                currentDice[nextToReroll] = RollDice(rng, sides);
                nextToReroll = (nextToReroll + 1) % numDice;
            }

            for (int d : currentDice)
                rollsRed[index] += d;
        }
    }

    // blue noise
    {
        // prime the dice
        std::vector<int> currentDice(numDice);
        for (int& d : currentDice)
            d = RollDice(rng, sides);

        int nextToReroll = 0;

        // do dice rolls
        std::vector<int> rolls(numRolls, 0);
        for (int index = 0; index < numRolls; ++index)
        {
            // reroll however many dice we should
            for (int rerollIndex = 0; rerollIndex < numRerolls; ++rerollIndex)
            {
                currentDice[nextToReroll] = RollDice(rng, sides);
                nextToReroll = (nextToReroll + 1) % numDice;
            }

            for (int i = 0; i < numDice; ++i)
                rollsBlue[index] += currentDice[i] * ((((index + i) % numDice) % 2) ? -1 : 1);
        }
    }

    std::vector<int> rolls(numRolls, 0);
    int minValue = 10000;
    int maxValue = -10000;
    for (int index = 0; index < numRolls; ++index)
    {
        rolls[index] = std::max(rollsRed[index], rollsBlue[index]);
        minValue = std::min(minValue, rolls[index]);
        maxValue = std::max(maxValue, rolls[index]);
    }

    // do tests
    char fileName[256];
    sprintf(fileName, "out/maxbluered%i", numDice);
    TestData(fileName, rolls, 6, minValue, maxValue);
}

int main(int argc, char** argv)
{
    int c_rolls[] =
    {
        16,
        128,
        1024,
        1024 * 16
    };

    for (int rollsIndex = 0; rollsIndex < _countof(c_rolls); ++rollsIndex)
    {
        int numRolls = c_rolls[rollsIndex];

        // white noise
        {
            // white noise drawn from rectangular distribution
            {
                RollDiceUncorrelated(6, numRolls);
            }

            // white noise drawn from a triangular distribution
            {
                RollDiceAdditiveUncorrelated(6, 2, numRolls);
                RollDiceAdditiveUncorrelated(6, 3, numRolls);
            }

            // white noise drawn from a triangular distribution
            {
                RollDiceSubtractUncorrelated(6, numRolls, 2);
                RollDiceSubtractUncorrelated(6, numRolls, 3);
            }

            // white noise drawn from a gaussian distribution
            {
                GaussianUncorrelated(0.0f, 1.0f, numRolls);
                GaussianUncorrelated(1.0f, 3.0f, numRolls);
            }
        }

        // red noise drawn from a triangular distribution
        {
            RollDiceAdditiveCorrelated(6, 2, numRolls, 1);
            RollDiceAdditiveCorrelated(6, 3, numRolls, 1);
            RollDiceAdditiveCorrelated(6, 4, numRolls, 1);
            RollDiceAdditiveCorrelated(6, 10, numRolls, 1);
            RollDiceAdditiveCorrelated(6, 20, numRolls, 1);
            RollDiceAdditiveCorrelated(6, 20, numRolls, 5);
            RollDiceAdditiveCorrelated(6, 20, numRolls, 10);
            RollDiceAdditiveCorrelated(6, 20, numRolls, 15);
            RollDiceAdditiveCorrelated(6, 20, numRolls, 19);
        }

        // red noise from a gaussian distribution
        GaussianAdditiveCorrelated(20, 0.0f, 1.0f, numRolls, 1);

        // blue noise drawn from a triangular distribution
        {
            RollDiceSubtractCorrelated(6, numRolls, 2, 1);
            RollDiceSubtractCorrelated(6, numRolls, 3, 1);
            RollDiceSubtractCorrelated(6, numRolls, 4, 1);
            RollDiceSubtractCorrelated(6, numRolls, 10, 1);
            RollDiceSubtractCorrelated(6, numRolls, 20, 1);
            RollDiceSubtractCorrelated(6, numRolls, 20, 5);
            RollDiceSubtractCorrelated(6, numRolls, 20, 10);
            RollDiceSubtractCorrelated(6, numRolls, 20, 15);
            RollDiceSubtractCorrelated(6, numRolls, 20, 19);
        }

        // blue noise from a gaussian distribution
        GaussianSubtractCorrelated(20, 0.0f, 1.0f, numRolls, 1);

        // add blue and red noise
        AddBlueRed(numRolls);

        // mix (lerp by a constant) red and blue noise
        MixBlueRed(numRolls, 0.25f);
        MixBlueRed(numRolls, 0.5f);
        MixBlueRed(numRolls, 0.75f);

        // max red and blue noise
        MaxBlueRed(numRolls);
    }

    return 0;
}


/*

NOTES:
* more dice rerolls = closer to white noise. So, fewer rerolls = more strongly filtered.
* more dice = better color and distribution
* more rolls = better color and distribution

Interesting info about how it relates to convolution
https://stats.stackexchange.com/questions/331973/why-is-the-sum-of-two-random-variables-a-convolution/331983#331983
https://twitter.com/ApoorvaJ/status/1144302805819305989?s=03

MaxBlueRed is a throwback to the maxing of uniform random values post

*/