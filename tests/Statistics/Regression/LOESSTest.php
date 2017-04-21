<?php
namespace MathPHP\Statistics\Regression;

class LOESSTest extends \PHPUnit_Framework_TestCase
{
    /**
     * @dataProvider dataProviderForEvaluate
     */
    public function testEvaluate(array $points, $α, $λ, $yHat)
    {
        $loess     = new LOESS($points, $α, $λ);
        $test_yHat = $loess->yHat();
        foreach ($test_yHat as $key => $value) {
            $this->assertEquals($yHat[$key], $value, '', .0000001);
        }
    }
    public function dataProviderForEvaluate()
    {
        
        // data from http://www.itl.nist.gov/div898/handbook/pmd/section1/dep/dep144.htm
        return [
            [
                [ [0.5578196, 18.63654], [2.0217271, 103.49646], [2.5773252, 150.35391], [3.4140288, 190.51031], [4.3014084, 208.70115], [4.7448394, 213.71135], [5.1073781, 228.49353], [6.5411662, 233.55387], [6.7216176, 234.55054], [7.2600583, 223.89225], [8.1335874, 227.68339], [9.1224379, 223.91982], [11.9296663, 168.01999], [12.3797674, 164.9575], [13.2728619, 152.61107], [14.2767453, 160.78742], [15.3731026, 168.55567], [15.6476637, 152.42658], [18.5605355, 221.70702], [18.5866354, 222.6904], [18.7572812, 243.18828] ],
                1/3,
                1,
                [20.5930233664173, 107.160307187003, 139.767381190249, 174.263043459836, 207.233382548945, 216.661586014446, 220.544479834053, 229.860693009198, 229.834712999942, 229.430115826698, 226.604459037782, 220.390409885032, 172.347999406849, 163.841661310125, 161.848970688604, 160.335083688228, 160.191989312566, 161.055592538928, 227.339955872967, 227.898534977536, 231.558556343728],
            ],
        ];
    }

    /**
     * @dataProvider dataProviderForSmoothnessParameterOutOfBoundsException
     */
    public function testSmoothnessParameterOutOfBoundsException(array $points, $α, $λ)
    {
        $this->setExpectedException('MathPHP\Exception\OutOfBoundsException');
        $loess = new LOESS($points, $α, $λ);
    }

    public function dataProviderForSmoothnessParameterOutOfBoundsException()
    {
        return [
            [
                [ [0.5578196, 18.63654], [2.0217271, 103.49646], [2.5773252, 150.35391], [3.4140288, 190.51031], [4.3014084, 208.70115], [4.7448394, 213.71135], [5.1073781, 228.49353], [6.5411662, 233.55387], [6.7216176, 234.55054], [7.2600583, 223.89225], [8.1335874, 227.68339], [9.1224379, 223.91982], [11.9296663, 168.01999], [12.3797674, 164.9575], [13.2728619, 152.61107], [14.2767453, 160.78742], [15.3731026, 168.55567], [15.6476637, 152.42658], [18.5605355, 221.70702], [18.5866354, 222.6904], [18.7572812, 243.18828] ],
                -50, 1
            ],
            [
                [ [0.5578196, 18.63654], [2.0217271, 103.49646], [2.5773252, 150.35391], [3.4140288, 190.51031], [4.3014084, 208.70115], [4.7448394, 213.71135], [5.1073781, 228.49353], [6.5411662, 233.55387], [6.7216176, 234.55054], [7.2600583, 223.89225], [8.1335874, 227.68339], [9.1224379, 223.91982], [11.9296663, 168.01999], [12.3797674, 164.9575], [13.2728619, 152.61107], [14.2767453, 160.78742], [15.3731026, 168.55567], [15.6476637, 152.42658], [18.5605355, 221.70702], [18.5866354, 222.6904], [18.7572812, 243.18828] ],
                1.1, 1
            ],
            [
                [ [0.5578196, 18.63654], [2.0217271, 103.49646], [2.5773252, 150.35391], [3.4140288, 190.51031], [4.3014084, 208.70115], [4.7448394, 213.71135], [5.1073781, 228.49353], [6.5411662, 233.55387], [6.7216176, 234.55054], [7.2600583, 223.89225], [8.1335874, 227.68339], [9.1224379, 223.91982], [11.9296663, 168.01999], [12.3797674, 164.9575], [13.2728619, 152.61107], [14.2767453, 160.78742], [15.3731026, 168.55567], [15.6476637, 152.42658], [18.5605355, 221.70702], [18.5866354, 222.6904], [18.7572812, 243.18828] ],
                3, 1
            ],
        ];
    }
}
