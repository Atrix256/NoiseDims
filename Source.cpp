
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

inline float Lerp(float a, float b, float t)
{
    return a * (1.0f - t) + b * t;
}

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

    // Zero out DC, we don't really care about it, and the value is huge.
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

void TestData(const char* baseFileName_, const std::vector<std::vector<float>>& tests, float mean, float sigma)
{
    char baseFileName[512];
    if (tests.size() > 1)
        sprintf(baseFileName, "%s_x%zu", baseFileName_, tests.size());
    else
        strcpy(baseFileName, baseFileName_);

    // get min / max values
    float minValue = tests[0][0];
    float maxValue = minValue;
    for (const std::vector<float>& test : tests)
    {
        for (float f : test)
        {
            minValue = std::min(minValue, f);
            maxValue = std::max(maxValue, f);
        }
    }

    // average the tests to make a histogram
    int numRolls = int(tests[0].size());
    std::vector<float> histogram(101, 0);
    {
        float normalizedHistogramCount = 1.0f / float(tests.size());
        for (int testIndex = 0; testIndex < tests.size(); ++testIndex)
        {
            for (int rollIndex = 0; rollIndex < numRolls; ++rollIndex)
            {
                size_t percent = size_t(100.0f * (tests[testIndex][rollIndex] - minValue) / (maxValue - minValue));
                histogram[percent] += normalizedHistogramCount;
            }
        }
    }

    // write histogram
    {
        char fileName[256];
        sprintf(fileName, "%s_%i_%i_%i_histogram.csv", baseFileName, int(mean), int(sigma), numRolls);
        FILE* file = nullptr;
        fopen_s(&file, fileName, "w+t");
        for (int index = 0; index < histogram.size(); ++index)
            fprintf(file, "\"%i\",\"%i\"\n", index, int(histogram[index]));
        fclose(file);
    }

    // for each test:
    // 1) DFT the rolls
    // 2) Adjust the average DFT
    std::vector<float> rollsFloatMagAvg(numRolls);
    {
        for (int testIndex = 0; testIndex < tests.size(); ++testIndex)
        {
            std::vector<float> rollsFloatMag;
            DFT(tests[testIndex], rollsFloatMag);

            // incremental average
            for (int rollIndex = 0; rollIndex < numRolls; ++rollIndex)
                rollsFloatMagAvg[rollIndex] = Lerp(rollsFloatMagAvg[rollIndex], rollsFloatMag[rollIndex], 1.0f / float(testIndex + 1));
        }

        char fileName[256];
        sprintf(fileName, "%s_%i_%i_%i_DFTMag.csv", baseFileName, int(mean), int(sigma), numRolls);
        FILE* file = nullptr;
        fopen_s(&file, fileName, "w+t");
        for (int index = 0; index < int(numRolls); ++index)
            fprintf(file, "\"%i\",\"%f\"\n", int(index) - int(numRolls) / 2, rollsFloatMagAvg[index]);
        fclose(file);
    }

    // write actual values from the first test
    {
        char fileName[256];
        sprintf(fileName, "%s_%i_%i_%i.txt", baseFileName, int(mean), int(sigma), numRolls);
        FILE* file = nullptr;
        fopen_s(&file, fileName, "w+t");
        for (float f : tests[0])
            fprintf(file, "%f\n", f);
        fclose(file);
    }
}

void TestData(const char* baseFileName_, const std::vector<std::vector<int>>& tests, int sides, int minValue, int maxValue)
{
    char baseFileName[512];
    if (tests.size() > 1)
        sprintf(baseFileName, "%s_x%zu", baseFileName_, tests.size());
    else
        strcpy(baseFileName, baseFileName_);

    // average the tests to make a histogram
    std::vector<float> histogram(maxValue - minValue + 1, 0);
    int numRolls = int(tests[0].size());
    {
        float normalizedHistogramCount = 1.0f / float(tests.size());
        for (int testIndex = 0; testIndex < tests.size(); ++testIndex)
        {
            for (int rollIndex = 0; rollIndex < numRolls; ++rollIndex)
                histogram[tests[testIndex][rollIndex] - minValue] += normalizedHistogramCount;
        }
    }

    // write histogram
    {
        char fileName[256];
        sprintf(fileName, "%s_%i_%i_histogram.csv", baseFileName, sides, numRolls);
        FILE* file = nullptr;
        fopen_s(&file, fileName, "w+t");
        for (int index = 0; index < histogram.size(); ++index)
            fprintf(file, "\"%i\",\"%i\"\n", index + minValue, int(histogram[index]));
        fclose(file);
    }

    // write DFT
    {
        // for each test:
        // 1) Make the rolls into floats
        // 2) DFT the rolls
        // 3) Adjust the average DFT
        std::vector<float> rollsFloatMagAvg(numRolls);
        {
            std::vector<float> rollsFloat(numRolls);
            for (int testIndex = 0; testIndex < tests.size(); ++testIndex)
            {
                for (int rollIndex = 0; rollIndex < numRolls; ++rollIndex)
                    rollsFloat[rollIndex] = float(tests[testIndex][rollIndex]);

                std::vector<float> rollsFloatMag;
                DFT(rollsFloat, rollsFloatMag);

                // incremental average
                for (int rollIndex = 0; rollIndex < numRolls; ++rollIndex)
                    rollsFloatMagAvg[rollIndex] = Lerp(rollsFloatMagAvg[rollIndex], rollsFloatMag[rollIndex], 1.0f / float(testIndex+1));
            }
        }

        char fileName[256];
        sprintf(fileName, "%s_%i_%i_DFTMag.csv", baseFileName, sides, numRolls);
        FILE* file = nullptr;
        fopen_s(&file, fileName, "w+t");
        for (int index = 0; index < int(numRolls); ++index)
            fprintf(file, "\"%i\",\"%f\"\n", int(index) - int(numRolls) / 2, rollsFloatMagAvg[index]);
        fclose(file);
    }

    // write actual values from the first test
    {
        char fileName[256];
        sprintf(fileName, "%s_%i_%i.txt", baseFileName, sides, numRolls);
        FILE* file = nullptr;
        fopen_s(&file, fileName, "w+t");
        for (int i : tests[0])
            fprintf(file, "%i\n", i);
        fclose(file);
    }
}

void RollDiceUncorrelated(int sides, int numRolls, int numTests)
{
    std::mt19937 rng(GetRNGSeed());

    // do each test
    std::vector<std::vector<int>> tests(numTests);
    for (int testIndex = 0; testIndex < numTests; ++testIndex)
    {
        // allocate space for this test
        std::vector<int>& rolls = tests[testIndex];
        rolls.resize(numRolls, 0);

        // do dice rolls
        for (int index = 0; index < numRolls; ++index)
            rolls[index] = RollDice(rng, sides);
    }

    // do tests
    TestData("out/uncorrelated", tests, sides, 0, 5);
}

void GaussianUncorrelated(float mean, float sigma, int numRolls, int numTests)
{
    std::mt19937 rng(GetRNGSeed());
    std::normal_distribution<float> dist(mean, sigma);

    // do each test
    std::vector<std::vector<float>> tests(numTests);
    for (int testIndex = 0; testIndex < numTests; ++testIndex)
    {
        // allocate space for this test
        std::vector<float>& rolls = tests[testIndex];
        rolls.resize(numRolls, 0);

        // do dice rolls
        for (int index = 0; index < numRolls; ++index)
            rolls[index] = dist(rng);
    }

    // do tests
    TestData("out/gaussianuncorrelated", tests, mean, sigma);
}

void RollDiceAdditiveUncorrelated(int sides, int numDice, int numRolls, int numTests)
{
    std::mt19937 rng(GetRNGSeed());

    // do each test
    std::vector<std::vector<int>> tests(numTests);
    for (int testIndex = 0; testIndex < numTests; ++testIndex)
    {
        // allocate space for this test
        std::vector<int>& rolls = tests[testIndex];
        rolls.resize(numRolls, 0);

        // do dice rolls
        for (int index = 0; index < numRolls; ++index)
        {
            for (int index2 = 0; index2 < numDice; ++index2)
                rolls[index] += RollDice(rng, sides);
        }
    }

    // do tests
    char fileName[256];
    sprintf(fileName, "out/adduncorrelated%i", numDice);
    TestData(fileName, tests, sides, 0, numDice * 5);
}

void RollDiceAdditiveCorrelated(int sides, int numDice, int numRolls, int numRerolls, int numTests)
{
    std::mt19937 rng(GetRNGSeed());

    // no funny business boyos
    numRerolls = std::min(numRerolls, numDice);

    // do each test
    std::vector<std::vector<int>> tests(numTests);
    for (int testIndex = 0; testIndex < numTests; ++testIndex)
    {
        // allocate space for this test
        std::vector<int>& rolls = tests[testIndex];
        rolls.resize(numRolls, 0);

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
                rolls[index] += d;
        }
    }

    // do tests
    char fileName[256];
    sprintf(fileName, "out/addcorrelated%i_%i", numDice, numRerolls);
    TestData(fileName, tests, sides, 0, numDice * 5);
}

void GaussianAdditiveCorrelated(int numDice, float mean, float sigma, int numRolls, int numRerolls, int numTests)
{
    std::mt19937 rng(GetRNGSeed());
    std::normal_distribution<float> dist(mean, sigma);

    // no funny business boyos
    numRerolls = std::min(numRerolls, numDice);

    // do each test
    std::vector<std::vector<float>> tests(numTests);
    for (int testIndex = 0; testIndex < numTests; ++testIndex)
    {
        // allocate space for this test
        std::vector<float>& rolls = tests[testIndex];
        rolls.resize(numRolls, 0);

        // prime the dice
        std::vector<float> currentDice(numDice);
        for (float& d : currentDice)
            d = dist(rng);

        int nextToReroll = 0;

        // do dice rolls
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
    }

    // do tests
    char fileName[256];
    sprintf(fileName, "out/gaussianaddcorrelated%i_%i", numDice, numRerolls);
    TestData(fileName, tests, mean, sigma);
}

void RollDiceSubtractUncorrelated(int sides, int numRolls, int numDice, int numTests)
{
    std::mt19937 rng(GetRNGSeed());

    // do each test
    std::vector<std::vector<int>> tests(numTests);
    for (int testIndex = 0; testIndex < numTests; ++testIndex)
    {
        // allocate space for this test
        std::vector<int>& rolls = tests[testIndex];
        rolls.resize(numRolls, 0);

        // do dice rolls
        std::vector<int> currentDice(numDice);
        for (int index = 0; index < numRolls; ++index)
        {
            for (int index2 = 0; index2 < numDice; ++index2)
                currentDice[index2] = RollDice(rng, sides);

            for (int i = 0; i < numDice; ++i)
                rolls[index] += currentDice[i] * ((((index + i) % numDice) % 2) ? -1 : 1);
        }
    }

    // do tests
    char fileName[256];
    sprintf(fileName, "out/subuncorrelated%i", numDice);
    TestData(fileName, tests, sides, -(sides - 1)*(numDice-1), (sides - 1)*(numDice-1));
}

void RollDiceSubtractCorrelated(int sides, int numRolls, int numDice, int numRerolls, int numTests)
{
    std::mt19937 rng(GetRNGSeed());

    // no funny business boyos
    numRerolls = std::min(numRerolls, numDice);

    // do each test
    std::vector<std::vector<int>> tests(numTests);
    for (int testIndex = 0; testIndex < numTests; ++testIndex)
    {
        // allocate space for this test
        std::vector<int>& rolls = tests[testIndex];
        rolls.resize(numRolls, 0);

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

            for (int i = 0; i < numDice; ++i)
                rolls[index] += currentDice[i] * ((((index + i) % numDice) % 2) ? -1 : 1);
        }
    }

    // do tests
    char fileName[256];
    sprintf(fileName, "out/subcorrelated%i_%i", numDice, numRerolls);
    TestData(fileName, tests, sides, -(sides - 1)*(numDice - 1), (sides - 1)*(numDice - 1));
}

void GaussianSubtractCorrelated(int numDice, float mean, float sigma, int numRolls, int numRerolls, int numTests)
{
    std::mt19937 rng(GetRNGSeed());
    std::normal_distribution<float> dist(mean, sigma);

    // no funny business boyos
    numRerolls = std::min(numRerolls, numDice);

    // do each test
    std::vector<std::vector<float>> tests(numTests);
    for (int testIndex = 0; testIndex < numTests; ++testIndex)
    {
        // allocate space for this test
        std::vector<float>& rolls = tests[testIndex];
        rolls.resize(numRolls, 0);

        // prime the dice
        std::vector<float> currentDice(numDice);
        for (float& d : currentDice)
            d = dist(rng);

        int nextToReroll = 0;

        // do dice rolls
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
    }

    // do tests
    char fileName[512];
    sprintf(fileName, "out/gaussiansubcorrelated%i_%i", numDice, numRerolls);
    TestData(fileName, tests, mean, sigma);
}

void AddBlueRed(int numRolls, int numTests)
{
    std::mt19937 rng(GetRNGSeed());

    int numRerolls = 1;
    int sides = 6;
    int numDice = 20;

    int minValue = 10000;
    int maxValue = -10000;

    // do each test
    std::vector<std::vector<int>> tests(numTests);
    for (int testIndex = 0; testIndex < numTests; ++testIndex)
    {
        // allocate space for this test
        std::vector<int>& rolls = tests[testIndex];
        rolls.resize(numRolls, 0);

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

        for (int index = 0; index < numRolls; ++index)
        {
            rolls[index] = rollsRed[index] + rollsBlue[index];
            minValue = std::min(minValue, rolls[index]);
            maxValue = std::max(maxValue, rolls[index]);
        }
    }

    // do tests
    char fileName[256];
    sprintf(fileName, "out/addbluered%i", numDice);
    TestData(fileName, tests, 6, minValue, maxValue);
}

void MixBlueRed(int numRolls, float redPercentage, int numTests)
{
    std::mt19937 rng(GetRNGSeed());
    float bluePercentage = 1.0f - redPercentage;

    int numRerolls = 1;
    int sides = 6;
    int numDice = 20;

    int minValue = 10000;
    int maxValue = -10000;

    // do each test
    std::vector<std::vector<int>> tests(numTests);
    for (int testIndex = 0; testIndex < numTests; ++testIndex)
    {
        // allocate space for this test
        std::vector<int>& rolls = tests[testIndex];
        rolls.resize(numRolls, 0);

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

        for (int index = 0; index < numRolls; ++index)
        {
            rolls[index] = int(float(rollsRed[index]) * redPercentage + float(rollsBlue[index]) * bluePercentage);
            minValue = std::min(minValue, rolls[index]);
            maxValue = std::max(maxValue, rolls[index]);
        }
    }

    // do tests
    char fileName[256];
    sprintf(fileName, "out/mixbluered%i_%i", numDice, int(100.0f * redPercentage));
    TestData(fileName, tests, 6, minValue, maxValue);
}

void MaxBlueRed(int numRolls, int numTests)
{
    std::mt19937 rng(GetRNGSeed());

    int numRerolls = 1;
    int sides = 6;
    int numDice = 20;

    int minValue = 10000;
    int maxValue = -10000;

    // do each test
    std::vector<std::vector<int>> tests(numTests);
    for (int testIndex = 0; testIndex < numTests; ++testIndex)
    {
        // allocate space for this test
        std::vector<int>& rolls = tests[testIndex];
        rolls.resize(numRolls, 0);

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

        for (int index = 0; index < numRolls; ++index)
        {
            rolls[index] = std::max(rollsRed[index], rollsBlue[index]);
            minValue = std::min(minValue, rolls[index]);
            maxValue = std::max(maxValue, rolls[index]);
        }
    }

    // do tests
    char fileName[256];
    sprintf(fileName, "out/maxbluered%i", numDice);
    TestData(fileName, tests, 6, minValue, maxValue);
}

void VossAlgorithm(int sides, int numDice, int numRolls, int numTests)
{
    std::mt19937 rng(GetRNGSeed());

    // do each test
    std::vector<std::vector<int>> tests(numTests);
    for (int testIndex = 0; testIndex < numTests; ++testIndex)
    {
        // allocate space for this test
        std::vector<int>& rolls = tests[testIndex];
        rolls.resize(numRolls, 0);

        // prime the dice
        std::vector<int> currentDice(numDice);
        for (int& d : currentDice)
            d = RollDice(rng, sides);

        // get the value of the first dice roll
        for (int d : currentDice)
            rolls[0] += d;

        // do dice rolls
        for (size_t rollIndex = 1; rollIndex < numRolls; ++rollIndex)
        {
            // re-roll the dice that need changing. A varying number each iteration.
            for (size_t diceIndex = 0; diceIndex < numDice; ++diceIndex)
            {
                size_t mask = size_t(1) << diceIndex;

                if ((rollIndex & mask) != ((rollIndex - 1) & mask))
                    currentDice[diceIndex] = RollDice(rng, sides);
            }

            // sum the dice
            for (int d : currentDice)
                rolls[rollIndex] += d;
        }
    }

    // do tests
    char fileName[256];
    sprintf(fileName, "out/voss%i_%i", numDice, sides);
    TestData(fileName, tests, sides, 0, numDice * (sides-1));
}

size_t LeadingZeros(size_t value)
{
    size_t ret = 0;
    while (value && (value & 1) == 0)
    {
        value /= 2;
        ret++;
    }
    return ret;
}

void VossMcCartneyAlgorithm(int sides, int numDice, int numRolls, int numTests)
{
    std::mt19937 rng(GetRNGSeed());

    // do each test
    std::vector<std::vector<int>> tests(numTests);
    for (int testIndex = 0; testIndex < numTests; ++testIndex)
    {
        // allocate space for this test
        std::vector<int>& rolls = tests[testIndex];
        rolls.resize(numRolls, 0);

        // prime the dice
        std::vector<int> currentDice(numDice);
        for (int& d : currentDice)
            d = RollDice(rng, sides);

        // get the value of the first dice roll
        for (int d : currentDice)
            rolls[0] += d;

        // do dice rolls
        for (size_t rollIndex = 1; rollIndex < numRolls; ++rollIndex)
        {
            // reroll one die each iteration. Constant workload.
            // NOTE: since diceToReroll can be out of bounds, i'm using modulus to keep it in bounds.  That doesn't feel correct but i'm not sure what a better option is.
            size_t diceToReroll = LeadingZeros(rollIndex) % numDice;
            currentDice[diceToReroll] = RollDice(rng, sides);

            // sum the dice
            for (int d : currentDice)
                rolls[rollIndex] += d;
        }
    }

    // do tests
    char fileName[256];
    sprintf(fileName, "out/vossmc%i_%i", numDice, sides);
    TestData(fileName, tests, sides, 0, numDice * (sides - 1));
}

void VossMcCartneyStochasticAlgorithm(int sides, int numDice, int numRolls, int numTests)
{
    std::mt19937 rng(GetRNGSeed());

    // do each test
    std::vector<std::vector<int>> tests(numTests);
    for (int testIndex = 0; testIndex < numTests; ++testIndex)
    {
        // allocate space for this test
        std::vector<int>& rolls = tests[testIndex];
        rolls.resize(numRolls, 0);

        // prime the dice
        std::vector<int> currentDice(numDice);
        for (int& d : currentDice)
            d = RollDice(rng, sides);

        // get the value of the first dice roll
        for (int d : currentDice)
            rolls[0] += d;

        // do dice rolls
        std::geometric_distribution<int> dist(0.5);
        for (size_t rollIndex = 1; rollIndex < numRolls; ++rollIndex)
        {
            // reroll one die each iteration. Constant workload.
            // NOTE: since diceToReroll can be out of bounds, i'm using modulus to keep it in bounds.  That doesn't feel correct but i'm not sure what a better option is.
            size_t diceToReroll = dist(rng) % numDice;
            currentDice[diceToReroll] = RollDice(rng, sides);

            // sum the dice
            for (int d : currentDice)
                rolls[rollIndex] += d;
        }
    }

    // do tests
    char fileName[256];
    sprintf(fileName, "out/vossmcstoch%i_%i", numDice, sides);
    TestData(fileName, tests, sides, 0, numDice * (sides - 1));
}

int main(int argc, char** argv)
{
    struct TestParams
    {
        int numRolls;
        int numTests;
    };

    TestParams c_testParams[] =
    {
        {    16,   1},
        {   128,   1},
        {  1024,   1},
        {1024*8,   1},
        {1024*8, 100},
    };

    for (int rollsIndex = 0; rollsIndex < _countof(c_testParams); ++rollsIndex)
    {
        int numRolls = c_testParams[rollsIndex].numRolls;
        int numTests = c_testParams[rollsIndex].numTests;

        // white noise
        {
            // white noise drawn from rectangular distribution
            {
                RollDiceUncorrelated(6, numRolls, numTests);
            }

            // white noise drawn from a triangular distribution
            {
                RollDiceAdditiveUncorrelated(6, 2, numRolls, numTests);
                RollDiceAdditiveUncorrelated(6, 3, numRolls, numTests);
                RollDiceAdditiveUncorrelated(6, 20, numRolls, numTests);
            }

            // white noise drawn from a triangular distribution
            {
                RollDiceSubtractUncorrelated(6, numRolls, 2, numTests);
                RollDiceSubtractUncorrelated(6, numRolls, 3, numTests);
                RollDiceSubtractUncorrelated(6, numRolls, 20, numTests);
            }

            // white noise drawn from a gaussian distribution
            {
                GaussianUncorrelated(0.0f, 1.0f, numRolls, numTests);
                GaussianUncorrelated(1.0f, 3.0f, numRolls, numTests);
            }
        }

        // red noise drawn from a triangular distribution
        {
            RollDiceAdditiveCorrelated(6, 2, numRolls, 1, numTests);
            RollDiceAdditiveCorrelated(6, 3, numRolls, 1, numTests);
            RollDiceAdditiveCorrelated(6, 4, numRolls, 1, numTests);
            RollDiceAdditiveCorrelated(6, 10, numRolls, 1, numTests);
            RollDiceAdditiveCorrelated(6, 20, numRolls, 1, numTests);
            RollDiceAdditiveCorrelated(6, 20, numRolls, 5, numTests);
            RollDiceAdditiveCorrelated(6, 20, numRolls, 10, numTests);
            RollDiceAdditiveCorrelated(6, 20, numRolls, 15, numTests);
            RollDiceAdditiveCorrelated(6, 20, numRolls, 19, numTests);
        }

        // red noise from a gaussian distribution
        GaussianAdditiveCorrelated(20, 0.0f, 1.0f, numRolls, 1, numTests);

        // blue noise drawn from a triangular distribution
        {
            RollDiceSubtractCorrelated(6, numRolls, 2, 1, numTests);
            RollDiceSubtractCorrelated(6, numRolls, 3, 1, numTests);
            RollDiceSubtractCorrelated(6, numRolls, 4, 1, numTests);
            RollDiceSubtractCorrelated(6, numRolls, 10, 1, numTests);
            RollDiceSubtractCorrelated(6, numRolls, 20, 1, numTests);
            RollDiceSubtractCorrelated(6, numRolls, 20, 5, numTests);
            RollDiceSubtractCorrelated(6, numRolls, 20, 10, numTests);
            RollDiceSubtractCorrelated(6, numRolls, 20, 15, numTests);
            RollDiceSubtractCorrelated(6, numRolls, 20, 19, numTests);
        }

        // blue noise from a gaussian distribution
        GaussianSubtractCorrelated(20, 0.0f, 1.0f, numRolls, 1, numTests);

        // add blue and red noise
        AddBlueRed(numRolls, numTests);

        // mix (lerp by a constant) red and blue noise
        MixBlueRed(numRolls, 0.01f, numTests);
        MixBlueRed(numRolls, 0.25f, numTests);
        MixBlueRed(numRolls, 0.5f, numTests);
        MixBlueRed(numRolls, 0.75f, numTests);
        MixBlueRed(numRolls, 0.99f, numTests);

        // max red and blue noise
        MaxBlueRed(numRolls, numTests);

        // pink noise
        {
            // pink noise with voss algorithm
            VossAlgorithm(6, 4, numRolls, numTests);

            // pink noise with McCartney modifications
            VossMcCartneyAlgorithm(6, 4, numRolls, numTests);

            // pink noise with stochastic modifications
            VossMcCartneyStochasticAlgorithm(6, 4, numRolls, numTests);
        }
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


Voss algorithm & another, well explained:
https://www.dsprelated.com/showarticle/908.php



MaxBlueRed is a throwback to the maxing of uniform random values post



NEXT:
* hyper gaussian distribution?
* box-muller transform blue noise?
* how to test how good it is... maybe fit a gaussian to the data points and compare the error vs the real gaussian that it's drawn from?

hyperuniformity stuff
https://en.wikipedia.org/wiki/Disordered_hyperuniformity
https://www.quantamagazine.org/hyperuniformity-found-in-birds-math-and-physics-20160712

correlated gaussian noise. but correlated how?
https://www.mathworks.com/matlabcentral/fileexchange/21156-correlated-gaussian-noise

*/