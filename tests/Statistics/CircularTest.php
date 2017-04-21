<?php
namespace MathPHP\Statistics;

class CircularTest extends \PHPUnit_Framework_TestCase
{
    /**
     * @dataProvider dataProviderForMean
     */
    public function testMean($angles, $mean)
    {
        $this->assertEquals($mean, Circular::mean($angles), '', 0.000001);
    }

    /**
     * Test data made with R package circular's function mean.circular()
     * https://cran.r-project.org/web/packages/circular/circular.pdf
     */
    public function dataProviderForMean()
    {
        $π = \M_PI;

        return [

            [[0, 2 * $π], 0],
            [[0, 0.5 * $π], 0.7853982],
            [[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], 0.5],
            [[0*$π, 0.1*$π, 0.2*$π, 0.3*$π, 0.4*$π, 0.5*$π, 0.6*$π, 0.7*$π, 0.8*$π, 0.9*$π, 1*$π], 1.570796],
            [[0, 0, 90], .5226276],
            [[1.4*$π, 1.7*$π, 1.75*$π, 2.54*$π, 4.32*$π], -0.4242655],
            [[5, 60, 340], -1.423654],
            [[5, 50, 150, 250], -0.9253517],
            [[10, 20, 30], -1.991149],
            [[355, 5, 15], -2.935443],

            // In this test case, we end up with
            // sin(0) + sin(π) = 0 + 0 = 0
            // cos(0) + cos(π) = 1 - 1 = 0
            // So it seems like it should end up as atan2(0, 0),
            // but since the sum of sins isn't perfectly 0, it is a very small floating point number,
            // like atan2(1.2246467991474E-16, 0),
            // which ends up as arctan(infinity) which equals 1.57079633.
            // R mean.circular results in NA,
            // but tested with Python scipi.stats.circmean(), it results in 1.5707963267948966,
            // which matches our PHP answer.
            [[0, $π], 1.5707963267948966],
        ];
    }

    /**
     * @dataProvider dataProviderForResultantLength
     */
    public function testResultantLength($angles, $mean)
    {
        $this->assertEquals($mean, Circular::resultantLength($angles), '', 0.00001);
    }

    /**
     * Test data made with custom R function:
     * resultantLength <- function(x) {
     *     sinSum = sum(sin(x))
     *     cosSum = sum(cos(x))
     *     R      = sqrt(sinSum^2 + cosSum^2)
     *     return(R)
     * }
     */
    public function dataProviderForResultantLength()
    {
        $π = \M_PI;

        return [
            [[0, $π], 1.224647e-16],
            [[0, 0.5, $π], 1],
            [[0, 2 * $π], 2],
            [[0, 0.5 * $π], 1.414214],
            [[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], 10.4581],
            [[0*$π, 0.1*$π, 0.2*$π, 0.3*$π, 0.4*$π, 0.5*$π, 0.6*$π, 0.7*$π, 0.8*$π, 0.9*$π, 1*$π], 6.313752],
            [[0, 0, 90], 1.791007],
            [[1.4*$π, 1.7*$π, 1.75*$π, 2.54*$π, 4.32*$π], 1.532213],
            [[5, 60, 340], 0.6201251],
            [[5, 50, 150, 250], 3.63869],
            [[10, 20, 30], 0.6781431],
            [[355, 5, 15], 1.507955],
        ];
    }

    /**
     * @dataProvider dataProviderForMeanResultantLength
     */
    public function testMeanResultantLength($angles, $mean)
    {
        $this->assertEquals($mean, Circular::meanResultantLength($angles), '', 0.000001);
    }

    /**
     * Test data made with custom R function:
     * meanResultantLength <- function(x) {
     *     n      = length(x)
     *     sinSum = sum(sin(x))
     *     cosSum = sum(cos(x))
     *     rho    = sqrt(sinSum^2 + cosSum^2) / n
     *     return(rho)
     * }
     */
    public function dataProviderForMeanResultantLength()
    {
        $π = \M_PI;

        return [
            [[0, $π], 6.123234e-17],
            [[0, 0.5, $π], 0.3333333],
            [[0, 2 * $π], 1],
            [[0, 0.5 * $π], 0.7071068],
            [[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], 0.9507365],
            [[0*$π, 0.1*$π, 0.2*$π, 0.3*$π, 0.4*$π, 0.5*$π, 0.6*$π, 0.7*$π, 0.8*$π, 0.9*$π, 1*$π], 0.5739774],
            [[0, 0, 90], 0.5970023],
            [[1.4*$π, 1.7*$π, 1.75*$π, 2.54*$π, 4.32*$π], 0.3064425],
            [[5, 60, 340], 0.2067084],
            [[5, 50, 150, 250], 0.9096725],
            [[10, 20, 30], 0.2260477],
            [[355, 5, 15], 0.5026515],
        ];
    }

    /**
     * @dataProvider dataProviderForVariance
     */
    public function testVariance($angles, $mean)
    {
        $this->assertEquals($mean, Circular::variance($angles), '', 0.000001);
    }

    /**
     * Test data made with R package circular's function var.circular()
     * https://cran.r-project.org/web/packages/circular/circular.pdf
     */
    public function dataProviderForVariance()
    {
        $π = \M_PI;

        return [
            [[0, $π], 1],
            [[0, 0.5, $π], 0.6666667],
            [[0, 2 * $π], 0],
            [[0, 0.5 * $π], 0.2928932],
            [[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], 0.04926349],
            [[0*$π, 0.1*$π, 0.2*$π, 0.3*$π, 0.4*$π, 0.5*$π, 0.6*$π, 0.7*$π, 0.8*$π, 0.9*$π, 1*$π], 0.4260226],
            [[0, 0, 90], 0.4029977],
            [[1.4*$π, 1.7*$π, 1.75*$π, 2.54*$π, 4.32*$π], 0.6935575],
            [[5, 60, 340], 0.7932916],
            [[5, 50, 150, 250], 0.09032747],
            [[10, 20, 30], 0.7739523],
            [[355, 5, 15], 0.4973485],
        ];
    }

    /**
     * @dataProvider dataProviderForStandardDeviation
     */
    public function testStandardDeviation($angles, $mean)
    {
        $this->assertEquals($mean, Circular::standardDeviation($angles), '', 0.000001);
    }

    /**
     * Test data made with R package circular's function sd.circular()
     * https://cran.r-project.org/web/packages/circular/circular.pdf
     */
    public function dataProviderForStandardDeviation()
    {
        $π = \M_PI;

        return [
            [[0, $π], 8.640817],
            [[0, 0.5, $π], 1.482304],
            [[0, 2 * $π], 0],
            [[0, 0.5 * $π], 0.8325546],
            [[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], 0.3178626],
            [[0*$π, 0.1*$π, 0.2*$π, 0.3*$π, 0.4*$π, 0.5*$π, 0.6*$π, 0.7*$π, 0.8*$π, 0.9*$π, 1*$π], 1.053722],
            [[0, 0, 90], 1.015711],
            [[1.4*$π, 1.7*$π, 1.75*$π, 2.54*$π, 4.32*$π], 1.538002],
            [[5, 60, 340], 1.775639],
            [[5, 50, 150, 250], 0.4351335],
            [[10, 20, 30], 1.724534],
            [[355, 5, 15], 1.172909],
        ];
    }

    public function testDescribe()
    {
        $stats = Circular::describe([5, 15, 355]);

        $this->assertTrue(is_array($stats));
        $this->assertArrayHasKey('n', $stats);
        $this->assertArrayHasKey('mean', $stats);
        $this->assertArrayHasKey('resultant_length', $stats);
        $this->assertArrayHasKey('mean_resultant_length', $stats);
        $this->assertArrayHasKey('variance', $stats);
        $this->assertArrayHasKey('sd', $stats);

        $this->assertTrue(is_int($stats['n']));
        $this->assertTrue(is_float($stats['mean']));
        $this->assertTrue(is_float($stats['resultant_length']));
        $this->assertTrue(is_float($stats['mean_resultant_length']));
        $this->assertTrue(is_float($stats['variance']));
        $this->assertTrue(is_float($stats['sd']));
    }
}
