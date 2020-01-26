
#include <catch2/catch.hpp>

#include "../src/stdafx.h"
#include "../src/solnp.hpp"

// Benchmark functions
#include "../src/benchmark/alkyla.hpp"
#include "../src/benchmark/box.hpp"
#include "../src/benchmark/entropy.hpp"
#include "../src/benchmark/powell.hpp"
#include "../src/benchmark/rosen_suzuki.hpp"
#include "../src/benchmark/wright_four.hpp"
#include "../src/benchmark/wright_nine.hpp"


TEST_CASE("Calculate the Alkyla function", "[alkyla]") {
    dlib::matrix<double, 10, 1> parameter_data;
    parameter_data =
            17.45,
            12.0,
            110.0,
            30.0,
            19.74,
            89.2,
            92.8,
            8.0,
            3.6,
            155.0;

    dlib::matrix<double, 8, 1> result = alkyla(parameter_data);

    CHECK(result(0) == Approx(-59.175999999999931));
    CHECK(result(1) == Approx(4.639999999999418));
    CHECK(result(2) == Approx(14.0));
    CHECK(result(3) == Approx(-58.999999999999773));
    CHECK(result(4) == Approx(1.015869200000000));
    CHECK(result(5) == Approx(0.997758620689655));
    CHECK(result(6) == Approx(0.391666666666666));
    CHECK(result(7) == Approx(0.938064516129032));


}


TEST_CASE("Optimize the Alkyla function manual hessian", "[alkyla]") {


    /* x0, lower, upper */
    dlib::matrix<double, 10, 3> parameter_data;
    parameter_data =
            17.45, 0.0, 20.0,
            12.0, 0.0, 16.0,
            110.0, 0.0, 120.0,
            30.0, 10.0, 50.0,
            19.74, 0.0, 20.0,
            89.2, 85.0, 93.0,
            92.8, 10.0, 95.0,
            8.0, 3.0, 12.0,
            3.6, 1.0, 4.0,
            155.0, 145.0, 162.0;

    /* Inequality function constraints.*/
    dlib::matrix<double, 4, 2> ib;
    ib =
            .99, 100.0 / 99.0,
            .99, 100.0 / 99.0,
            .9, 10.0 / 9.0,
            .99, 100.0 / 99.0;

    dlib::matrix<double> hessian_matrix;
    hessian_matrix = dlib::identity_matrix<double>(parameter_data.nr() + ib.nr());

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(alkyla_functor(), parameter_data, ib, hessian_matrix, logger, 0.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Validate values
    CHECK(result(0) == Approx(0.169963765876488e2));
    CHECK(result(1) == Approx(0.159994026791621e2));
    CHECK(result(2) == Approx(0.576883584245853e2));
    CHECK(result(3) == Approx(0.303248903549694e2));
    CHECK(result(4) == Approx(0.199999896454138e2));
    CHECK(result(5) == Approx(0.905654248087076e2));
    CHECK(result(6) == Approx(0.949999927142590e2));
    CHECK(result(7) == Approx(0.105901405233357e2));
    CHECK(result(8) == Approx(0.015616462840774e2));
    CHECK(result(9) == Approx(1.535353201975077e2));

    REQUIRE(calculate <= Approx(-1.726412486481025e2));

}


TEST_CASE("Calculate the Box function", "[box]") {
    dlib::matrix<double, 3, 1> parameter_data;
    parameter_data =
            1.1,
            1.1,
            9.0;

    dlib::matrix<double, 2, 1> result = box(parameter_data);

    CHECK(result(0) == Approx(-10.890000000000002));
    CHECK(result(1) == Approx(-55.560000000000002));
}

TEST_CASE("Optimize the Box function (case a)", "[box]") {


    /* x0 */
    dlib::matrix<double, 3, 3> parameter_data;
    parameter_data =
            1.1, 1.0, 10.0,
            1.1, 1.0, 10.0,
            9.0, 1.0, 10.0;

    dlib::matrix<double, 0, 0> ib;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(box_functor(), parameter_data, ib, logger, 1.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Check the parameters
    CHECK(result(0) == Approx(2.886775069536727));
    CHECK(result(1) == Approx(2.886775072009683));
    CHECK(result(2) == Approx(5.773407750048355));

    REQUIRE(calculate <= -48.112522068150462);

}


TEST_CASE("Optimize the Box function (case b)", "[box]") {

    /* x0 */
    dlib::matrix<double, 3, 3> parameter_data;
    parameter_data =
            5.5, 1.0, 10.0,
            5.5, 1.0, 10.0,
            5.5, 1.0, 10.0;

    dlib::matrix<double, 0, 0> ib;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(box_functor(), parameter_data, ib, logger, 1.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Check the parameters
    CHECK(result(0) == Approx(2.888765743268910));
    CHECK(result(1) == Approx(2.888765747765645));
    CHECK(result(2) == Approx(5.765448483893261));

    REQUIRE(calculate <= Approx(-48.112480408240664));

}


TEST_CASE("Calculate the Entropy function", "[entropy]") {
    dlib::matrix<double, 10, 1> parameter_data;
    parameter_data = 0.8474, 0.4524, 0.8075, 0.4832, 0.6135, 0.2749, 0.8807, 0.6538, 0.4899, 0.7741;

    dlib::matrix<double, 2, 1> result = entropy(parameter_data);

    CHECK(result(0) == Approx(4.849345605));
    CHECK(result(1) == Approx(-3.7226));

}


TEST_CASE("Optimize the Entropy function", "[entropy]") {


    /* x0, lower, upper */
    dlib::matrix<double, 10, 3> parameter_data;
    parameter_data =
            0.8474, 0.0, 10.0,
            0.4524, 0.0, 10.0,
            0.8075, 0.0, 10.0,
            0.4832, 0.0, 10.0,
            0.6135, 0.0, 10.0,
            0.2749, 0.0, 10.0,
            0.8807, 0.0, 10.0,
            0.6538, 0.0, 10.0,
            0.4899, 0.0, 10.0,
            0.7741, 0.0, 10.0;

    dlib::matrix<double, 0, 0> ib;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(entropy_functor(), parameter_data, ib, logger, 1.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Validate values
    CHECK(result(0) == Approx(0.857717373389226));
    CHECK(result(1) == Approx(0.857627630428985));
    CHECK(result(2) == Approx(0.858758912629577));
    CHECK(result(3) == Approx(0.857611377467009));
    CHECK(result(4) == Approx(0.857445851804114));
    CHECK(result(5) == Approx(0.857975885723616));
    CHECK(result(6) == Approx(2.279903520424556));
    CHECK(result(7) == Approx(0.857328700404536));
    CHECK(result(8) == Approx(0.857608122266489));
    CHECK(result(9) == Approx(0.858022625546437));

    REQUIRE(calculate <= Approx(0.185478885901993));

}


TEST_CASE("Calculate the Powell function", "[powell]") {
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data =
            -2.0,
            2.0,
            2.0,
            -1,
            -1;

    dlib::matrix<double, 4, 1> result = powell(parameter_data);

    CHECK(result(0) == Approx(0.000335462627903));
    CHECK(result(1) == Approx(4.0));
    CHECK(result(2) == Approx(-1.0));
    CHECK(result(3) == Approx(1.0));


}

TEST_CASE("Optimize the Powell function (rho == 0)", "[powell]") {


    /* x0 */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data =
            -2.0,
            2.0,
            2.0,
            -1,
            -1;

    dlib::matrix<double, 0, 0> ib;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(powell_functor(), parameter_data, ib, logger, 0.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Check the parameters
    CHECK(result(0) == Approx(-1.717126203723513));
    CHECK(result(1) == Approx(1.595689596511580));
    CHECK(result(2) == Approx(1.827278075550860));
    CHECK(result(3) == Approx(-0.763645042210886));
    CHECK(result(4) == Approx(-0.763645042234952));

    REQUIRE(calculate <= 0.053949827793391);

}

TEST_CASE("Optimize the Powell function (rho == 1)", "[powell]") {


    /* x0 */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data =
            -2.0,
            2.0,
            2.0,
            -1,
            -1;

    dlib::matrix<double, 0, 0> ib;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(powell_functor(), parameter_data, ib, logger, 1.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Check the parameters
    CHECK(result(0) == Approx(-1.717142506313303));
    CHECK(result(1) == Approx(1.595708459713134));
    CHECK(result(2) == Approx(1.827247731350245));
    CHECK(result(3) == Approx(-0.763643197991088));
    CHECK(result(4) == Approx(-0.763643197980140));

    REQUIRE(calculate <= 0.053949846871732);

}


TEST_CASE("Calculate the Wright4 function", "[wright4]") {
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = -1.0, 3.0, -0.5, -2.0, -3.0;

    dlib::matrix<double, 4, 1> result = wright_four(parameter_data);


    CHECK(result(0) == Approx(68.9375));
    CHECK(result(1) == Approx(1.632359313));
    CHECK(result(2) == Approx(-0.07842712475));
    CHECK(result(3) == Approx(1.0));

}


TEST_CASE("Optimize the Wright4 function (case a, rho==10)", "[wright4]") {


    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = 1, 1, 1, 1, 1;

    dlib::matrix<double, 0, 0> ib;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 10.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Validate values
    CHECK(result(0) == Approx(1.116639595144975));
    CHECK(result(1) == Approx(1.220448420527845));
    CHECK(result(2) == Approx(1.537782093973876));
    CHECK(result(3) == Approx(1.972752470314671));
    CHECK(result(4) == Approx(1.791088179957703));

    REQUIRE(calculate <= Approx(0.029310831271171));

}


TEST_CASE("Optimize the Wright4 function (case a, rho==1)", "[wright4]") {


    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = 1, 1, 1, 1, 1;

    dlib::matrix<double, 0, 0> ib;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 1.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Validate values
    CHECK(result(0) == Approx(1.116643803185402));
    CHECK(result(1) == Approx(1.220453096801827));
    CHECK(result(2) == Approx(1.537779890421572));
    CHECK(result(3) == Approx(1.972741018934506));
    CHECK(result(4) == Approx(1.791081456246007));

    REQUIRE(calculate <= Approx(0.029310831002758));

}


TEST_CASE("Optimize the Wright4 function (case a, rho==0)", "[wright4]") {


    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = 1, 1, 1, 1, 1;

    dlib::matrix<double, 0, 0> ib;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 0.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Validate values
    CHECK(result(0) == Approx(1.116609971954848));
    CHECK(result(1) == Approx(1.220440489315363));
    CHECK(result(2) == Approx(1.537788998574534));
    CHECK(result(3) == Approx(1.972781636643539));
    CHECK(result(4) == Approx(1.791135716667731));

    REQUIRE(calculate <= Approx(0.029310831942731));

}

TEST_CASE("Optimize the Wright4 function (case b, rho==10)", "[wright4]") {


    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = 2, 2, 2, 2, 2;

    dlib::matrix<double, 0, 0> ib;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 10.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Validate values
    CHECK(result(0) == Approx(1.116634078024861));
    CHECK(result(1) == Approx(1.220444234062229));
    CHECK(result(2) == Approx(1.537784320405761));
    CHECK(result(3) == Approx(1.972763493598540));
    CHECK(result(4) == Approx(1.791097035073237));

    REQUIRE(calculate <= Approx(0.029310831022048));

}


TEST_CASE("Optimize the Wright4 function (case b, rho==1)", "[wright4]") {


    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = 2, 2, 2, 2, 2;

    dlib::matrix<double, 0, 0> ib;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 1.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Validate values
    CHECK(result(0) == Approx(1.116634870722612));
    CHECK(result(1) == Approx(1.220440576122608));
    CHECK(result(2) == Approx(1.537785458091246));
    CHECK(result(3) == Approx(1.972770662225987));
    CHECK(result(4) == Approx(1.791095775742904));

    REQUIRE(calculate <= Approx(0.029310830648204));

}


TEST_CASE("Optimize the Wright4 function (case b, rho==0)", "[wright4]") {


    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = 2, 2, 2, 2, 2;

    dlib::matrix<double, 0, 0> ib;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 0.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Validate values
    CHECK(result(0) == Approx(1.116636615308806));
    CHECK(result(1) == Approx(1.220440910051733));
    CHECK(result(2) == Approx(1.537785097114965));
    CHECK(result(3) == Approx(1.972769218398222));
    CHECK(result(4) == Approx(1.791092977407627));

    REQUIRE(calculate <= Approx(0.029310830686406));

}

TEST_CASE("Optimize the Wright4 function (case c, rho==10)", "[wright4]") {


    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = -1, 3, -0.5, -2, -3;

    dlib::matrix<double, 0, 0> ib;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 10.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Validate values
    CHECK(result(0) == Approx(-0.703068933803915));
    CHECK(result(1) == Approx(2.635653020521741));
    CHECK(result(2) == Approx(-0.099089129462807));
    CHECK(result(3) == Approx(-1.797457648464959));
    CHECK(result(4) == Approx(-2.844671274183172));

    REQUIRE(calculate <= Approx(44.022877145171257));

}


TEST_CASE("Optimize the Wright4 function (case c, rho==1)", "[wright4]") {


    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = -1, 3, -0.5, -2, -3;

    dlib::matrix<double, 0, 0> ib;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 1.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Validate values
    CHECK(result(0) == Approx(-0.703523524834065));
    CHECK(result(1) == Approx(2.635727880354312));
    CHECK(result(2) == Approx(-0.096451726879527));
    CHECK(result(3) == Approx(-1.797997954636873));
    CHECK(result(4) == Approx(-2.842832966138447));

    REQUIRE(calculate <= Approx(44.022075061138295));

}


TEST_CASE("Optimize the Wright4 function (case c, rho==0)", "[wright4]") {


    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = -1, 3, -0.5, -2, -3;

    dlib::matrix<double, 0, 0> ib;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 0.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Validate values
    CHECK(result(0) == Approx(-0.703376376169597));
    CHECK(result(1) == Approx(2.635700821034921));
    CHECK(result(2) == Approx(-0.096657075140750));
    CHECK(result(3) == Approx(-1.797935433028946));
    CHECK(result(4) == Approx(-2.843427812167240));

    REQUIRE(calculate <= Approx(44.022128023467303));

}


TEST_CASE("Optimize the Wright4 function (case d, rho==10)", "[wright4]") {


    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = -1, 2, 1, -2, -2;

    dlib::matrix<double, 0, 0> ib;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 10, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Validate values
    CHECK(result(0) == Approx(-1.273052029422237));
    CHECK(result(1) == Approx(2.410353869751445));
    CHECK(result(2) == Approx(1.194859244564641));
    CHECK(result(3) == Approx(-0.154238130707192));
    CHECK(result(4) == Approx(-1.571027698602370));

    REQUIRE(calculate <= Approx(27.871905223431018));

}


TEST_CASE("Optimize the Wright4 function (case d, rho==1)", "[wright4]") {


    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = -1, 2, 1, -2, -2;

    dlib::matrix<double, 0, 0> ib;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 1.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Validate values
    CHECK(result(0) == Approx(-1.272710813834789));
    CHECK(result(1) == Approx(2.410353128557497));
    CHECK(result(2) == Approx(1.194780438668736));
    CHECK(result(3) == Approx(-0.154425728750688));
    CHECK(result(4) == Approx(-1.571448524460437));

    REQUIRE(calculate <= Approx(27.871903584038883));

}


TEST_CASE("Optimize the Wright4 function (case d, rho==0)", "[wright4]") {


    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = -1, 2, 1, -2, -2;

    dlib::matrix<double, 0, 0> ib;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 0.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Validate values
    CHECK(result(0) == Approx(-1.273023411221166));
    CHECK(result(1) == Approx(2.410362154910699));
    CHECK(result(2) == Approx(1.194843255003828));
    CHECK(result(3) == Approx(-0.154284646055543));
    CHECK(result(4) == Approx(-1.571062971404984));

    REQUIRE(calculate <= Approx(27.871904866028800));

}


TEST_CASE("Calculate the Wright9 function", "[wright9]") {
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data =
            1.091,
            -3.174,
            1.214,
            -1.614,
            2.134;

    dlib::matrix<double, 4, 1> result = wright_nine(parameter_data);

    CHECK(result(0) == Approx(-1.815401107e+3));
    CHECK(result(1) == Approx(0.019897305e+3));
    CHECK(result(2) == Approx(-0.001999274866e+3));
    CHECK(result(3) == Approx(0.007022058536e+3));


}


TEST_CASE("Optimize the Wright9 function (case a, rho==1)", "[wright9]") {


    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data =
            1,
            1,
            1,
            1,
            1;

    /* Inequality function constraints.*/
    dlib::matrix<double, 3, 2> ib;
    ib =
            -100, 20,
            -2, 100,
            5, 100;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(wright_nine_functor(), parameter_data, ib, logger, 1.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Validate values
    CHECK(result(0) == Approx(-0.082010879422399));
    CHECK(result(1) == Approx(3.692422439415791));
    CHECK(result(2) == Approx(2.487343192052682));
    CHECK(result(3) == Approx(0.377176847262690));
    CHECK(result(4) == Approx(0.173650632156765));

    REQUIRE(calculate <= Approx(-2.104078394423900e2));

}


TEST_CASE("Optimize the Wright9 function (case a, rho==100)", "[wright9]") {


    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data =
            1,
            1,
            1,
            1,
            1;

    /* Inequality function constraints.*/
    dlib::matrix<double, 3, 2> ib;
    ib =
            -100, 20,
            -2, 100,
            5, 100;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(wright_nine_functor(), parameter_data, ib, logger, 100.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Validate values
    CHECK(result(0) == Approx(-0.081246392868120));
    CHECK(result(1) == Approx(3.689847468079701));
    CHECK(result(2) == Approx(2.491144012826899));
    CHECK(result(3) == Approx(0.377589784583978));
    CHECK(result(4) == Approx(0.173399403653944));

    REQUIRE(calculate <= Approx(-2.104073432066561e2));

}


TEST_CASE("Optimize the Wright9 function (case b, rho==100)", "[wright9]") {


    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data =
            1.091,
            -3.174,
            1.214,
            -1.614,
            2.134;

    /* Inequality function constraints.*/
    dlib::matrix<double, 3, 2> ib;
    ib =
            -100, 20,
            -2, 100,
            5, 100;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(wright_nine_functor(), parameter_data, ib, logger, 100.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Validate values
    CHECK(result(0) == Approx(1.479634676414054));
    CHECK(result(1) == Approx(-2.636607129671721));
    CHECK(result(2) == Approx(1.054666762278611));
    CHECK(result(3) == Approx(-1.611508943269783));
    CHECK(result(4) == Approx(2.673892424752704));

    REQUIRE(calculate <= Approx(-2.500584227790517e3));

}


TEST_CASE("Calculate the Rosen-Suzuki function", "[rosen_suzuki]") {
    dlib::matrix<double, 4, 1> parameter_data;
    parameter_data = 1, 1, 1, 1;

    dlib::matrix<double, 4, 1> result = rosen_suzuki(parameter_data);

    CHECK(result(0) == Approx(-19.0));
    CHECK(result(1) == Approx(4.0));
    CHECK(result(2) == Approx(6.0));
    CHECK(result(3) == Approx(1.0));


}

TEST_CASE("Optimize the Rosen-Suzuki function", "[rosen_suzuki]") {


    /* x0 */
    dlib::matrix<double, 4, 1> parameter_data;
    parameter_data = 1, 1, 1, 1;

    /* Inequality function constraints.*/
    dlib::matrix<double, 3, 2> ib;
    ib =
            0, 1000,
            0, 1000,
            0, 1000;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(rosen_suzuki_functor(), parameter_data, ib, logger, 1.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Validate values
    CHECK(result(0) == Approx(0.000230374253029));
    CHECK(result(1) == Approx(0.998564179709879));
    CHECK(result(2) == Approx(2.000278419459943));
    CHECK(result(3) == Approx(-0.999859532645383));

    REQUIRE(calculate <= Approx(-43.999759237182886));

}