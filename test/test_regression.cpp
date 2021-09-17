/**
 * StreamHist testing module.
 * Most tests in this module are ported/adapted from the Clojure tests developed
 * for BigMl's "Streaming Histograms for Clojure/Java" [1].
 *
 * References
 * ----------
 * [1] https://github.com/bigmlcom/histogram
 *
 *
 * Copyright © 2015 Carson Farmer <carsonfarmer@gmail.com>
 * Copyright © 2013, 2014, 2015 BigML
 * Licensed under the Apache License, Version 2.0
 */


#include "streamhist/streamhist.h"
#include <catch2/catch.hpp>

#include <random>


using namespace streamhist;
using namespace std;


using Hist = streamhist::StreamHist<double>;


/*inline static std::vector<double> make_normal(size_t size) noexcept {
    std::random_device         rd;         // Will be used to obtain a seed for the random number engine
    std::mt19937               gen(rd());  // Standard mersenne_twister_engine seeded with rd()
    std::normal_distribution<> dis(0., 1.);

    std::vector<double> res;
    res.resize(size);

    for (size_t i = 0; i < size; i++) {
        res.push_back(dis(gen));
    }

    return res;
}*/

inline static bool about(double v1, double v2, double epsilon) noexcept {
    return std::abs(v1 - v2) <= epsilon;
}

TEST_CASE("example_part1", "StreamhistRegression") {

    /*randomcpp::seed(1700);
    auto data = make_normal(10000);*/

    // based on examples.ipynb
    // ensure the exact same values as pythons random values

    /**
     * normal_data = np.random.normal(0.0, 1.0, 10000)
     * print(normal_data)
     * with open('c:/temp/data.txt', 'w') as f:
     *     f.write("{\n")
     *     for item in normal_data:
     *         f.write("%s,\n" % item)
     *     f.write("};\n")
     */


#include "test_data3.h"

    {
        // h1 = StreamHist(maxbins=32)  # Create histogram with 32 bins
        // h1.update(normal_data)       # Add points all at once

        Hist h1(32);
        h1.update(data);

        REQUIRE(h1.sum(0.0) == Approx(5013.251966133985));
        REQUIRE(h1.density(0.0) == Approx(4029.1834532721346));

        auto d = h1.describe();
        /*for (auto& p : d) {
            std::cout << p.first << ": " << p.second << "\n";
        }*/

        REQUIRE(about(d["min"], -3.547153940496022, 0.05));
        REQUIRE(about(d["max"], 3.944891609124639, 0.05));
        REQUIRE(about(d["mean"], 0.00021621187019188817, 0.05));
        REQUIRE(about(d["var"], 0.9776676445197616, 0.05));
        REQUIRE(about(d["25%"], -0.6716281267050649, 0.05));
        REQUIRE(about(d["50%"], -0.0031658326488016794, 0.05));
        REQUIRE(about(d["75%"], 0.662280691581657, 0.05));
        REQUIRE(d["count"] == sizeof(data) / sizeof(data[0]));


        // h1.quantiles(0.5, 0.95, 0.99);  # Supports multiple quantile inputs
        double expected[] = { -0.0031658326488016794, 1.6623552691815082, 2.364513352627507 };
        size_t i          = 0;
        for (auto q : h1.quantiles(0.5, 0.95, 0.99)) {
            REQUIRE(about(q, expected[i++], 0.05));
        }
    }


    {
        // h1 = StreamHist(maxbins=32, weighted=True)  # Create histogram with 32 bins
        // h1.update(normal_data)       # Add points all at once

        Hist h1(32, true);
        h1.update(data);

        REQUIRE(h1.sum(0.0) == Approx(5018.37740081585));
        REQUIRE(h1.density(0.0) == Approx(4085.420804192544));

        auto d = h1.describe();
        /*for (auto& p : d) {
            std::cout << p.first << ": " << p.second << "\n";
        }*/

        REQUIRE(about(d["min"], -3.547153940496022, 0.05));
        REQUIRE(about(d["max"], 3.944891609124639, 0.05));
        REQUIRE(about(d["mean"], 0.00021621187019163984, 0.05));
        REQUIRE(about(d["var"], 0.9791353366986024, 0.05));
        REQUIRE(about(d["25%"], -0.6719326242816118, 0.05));
        REQUIRE(about(d["50%"], -0.004374895469705084, 0.05));
        REQUIRE(about(d["75%"], 0.6620785849105835, 0.05));
        REQUIRE(d["count"] == sizeof(data) / sizeof(data[0]));


        // h1.quantiles(0.5, 0.95, 0.99);  # Supports multiple quantile inputs
        double expected[] = { -0.004374895469705084, 1.6624420718883346, 2.3674215135178627 };
        size_t i          = 0;
        for (auto q : h1.quantiles(0.5, 0.95, 0.99)) {
            REQUIRE(about(q, expected[i++], 0.05));
        }
    }
}


TEST_CASE("example_part2", "StreamhistRegression") {

    /*randomcpp::seed(1700);
    auto data = make_normal(10000);*/

    // based on examples.ipynb
    // ensure the exact same values as pythons random values

    /**
     * mixed_normal_data = np.concatenate((
     *  np.random.normal(0.0, 0.2, 16000),
     *  np.random.normal(1.0, 0.2, 8000),
     *  np.random.normal(2.0, 0.2, 4000),
     *  np.random.normal(3.0, 0.2, 2000)
     *  ))
     * np.random.shuffle(mixed_normal_data)
     *
     * print(mixed_normal_data)
     * with open('c:/temp/data.txt', 'w') as f:
     *     f.write("{\n")
     *     for item in mixed_normal_data:
     *         f.write("%s,\n" % item)
     *     f.write("};\n")
     */

#include "test_data4.h"

    {
        // h3 = StreamHist(maxbins=8).update(mixed_normal_data)
        //  h4 = StreamHist(maxbins=64).update(mixed_normal_data)

        Hist h3(8);
        h3.update(data);
        Hist h4(32);
        h4.update(data);

        {
            auto d = h3.describe();
            /*for (auto& p : d) {
                std::cout << p.first << ": " << p.second << "\n";
            }*/

            REQUIRE(about(d["min"], -0.9224294154651361, 0.05));
            REQUIRE(about(d["max"], 3.7228156332509026, 0.05));
            REQUIRE(about(d["mean"], 0.7332564780333407, 0.05));
            REQUIRE(about(d["var"], 0.8717589972831745, 0.05));
            REQUIRE(about(d["25%"], -0.025022752392304337, 0.05));
            REQUIRE(about(d["50%"], 0.4181022716231957, 0.05));
            REQUIRE(about(d["75%"], 1.2255418397607567, 0.05));
            REQUIRE(d["count"] == sizeof(data) / sizeof(data[0]));
        }

        {
            auto d = h4.describe();
            /*for (auto& p : d) {
                std::cout << p.first << ": " << p.second << "\n";
            }*/

            REQUIRE(about(d["min"], -0.9224294154651361, 0.05));
            REQUIRE(about(d["max"], 3.7228156332509026, 0.05));
            REQUIRE(about(d["mean"], 0.7332564780333409, 0.05));
            REQUIRE(about(d["var"], 0.9035000864442686, 0.05));
            REQUIRE(about(d["25%"], -0.013509716909020641, 0.05));
            REQUIRE(about(d["50%"], 0.3065276859399, 0.05));
            REQUIRE(about(d["75%"], 1.1761322478033631, 0.05));
            REQUIRE(d["count"] == sizeof(data) / sizeof(data[0]));
        }
    }
}


TEST_CASE("iris_regression", "StreamhistRegression") {
    double sepal_length[]
        = { 5.1, 4.9, 4.7, 4.6, 5.0, 5.4, 4.6, 5.0, 4.4, 4.9, 5.4, 4.8, 4.8, 4.3, 5.8, 5.7, 5.4, 5.1, 5.7, 5.1, 5.4, 5.1, 4.6, 5.1, 4.8,
            5.0, 5.0, 5.2, 5.2, 4.7, 4.8, 5.4, 5.2, 5.5, 4.9, 5.0, 5.5, 4.9, 4.4, 5.1, 5.0, 4.5, 4.4, 5.0, 5.1, 4.8, 5.1, 4.6, 5.3, 5.0,
            7.0, 6.4, 6.9, 5.5, 6.5, 5.7, 6.3, 4.9, 6.6, 5.2, 5.0, 5.9, 6.0, 6.1, 5.6, 6.7, 5.6, 5.8, 6.2, 5.6, 5.9, 6.1, 6.3, 6.1, 6.4,
            6.6, 6.8, 6.7, 6.0, 5.7, 5.5, 5.5, 5.8, 6.0, 5.4, 6.0, 6.7, 6.3, 5.6, 5.5, 5.5, 6.1, 5.8, 5.0, 5.6, 5.7, 5.7, 6.2, 5.1, 5.7,
            6.3, 5.8, 7.1, 6.3, 6.5, 7.6, 4.9, 7.3, 6.7, 7.2, 6.5, 6.4, 6.8, 5.7, 5.8, 6.4, 6.5, 7.7, 7.7, 6.0, 6.9, 5.6, 7.7, 6.3, 6.7,
            7.2, 6.2, 6.1, 6.4, 7.2, 7.4, 7.9, 6.4, 6.3, 6.1, 7.7, 6.3, 6.4, 6.0, 6.9, 6.7, 6.9, 5.8, 6.8, 6.7, 6.7, 6.3, 6.5, 6.2, 5.9 };

    Hist h(32);
    h.update(sepal_length);

    Hist::Bin reg[] = { { 4.3, 1 },
                        { 4.425000000000001, 4 },
                        { 4.6, 4 },
                        { 4.771428571428571, 7 },
                        { 4.9625, 16 },
                        { 5.1, 9 },
                        { 5.2, 4 },
                        { 5.3, 1 },
                        { 5.3999999999999995, 6 },
                        { 5.5, 7 },
                        { 5.6000000000000005, 6 },
                        { 5.7, 8 },
                        { 5.8, 7 },
                        { 5.900000000000001, 3 },
                        { 6.0, 6 },
                        { 6.1000000000000005, 6 },
                        { 6.2, 4 },
                        { 6.299999999999999, 9 },
                        { 6.3999999999999995, 7 },
                        { 6.5, 5 },
                        { 6.6, 2 },
                        { 6.700000000000001, 8 },
                        { 6.8, 3 },
                        { 6.9, 4 },
                        { 7.0, 1 },
                        { 7.1, 1 },
                        { 7.2, 3 },
                        { 7.3, 1 },
                        { 7.4, 1 },
                        { 7.6, 1 },
                        { 7.7, 4 },
                        { 7.9, 1 } };

    REQUIRE(h.bins == reg);
}
