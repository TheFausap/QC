<?php
/**
 * Created by PhpStorm.
 * User: vandine
 * Date: 20/04/2017
 * Time: 20:51
 */
require_once(__DIR__ . '\..\vendor\autoload.php');

use MathPHP\LinearAlgebra\Matrix;
use MathPHP\LinearAlgebra\Vector;
use MathPHP\Number\Complex;

function tensorProduct(Vector $v1, Vector $v2)
{
    $lv1 = $v1->getN();
    $lv2 = $v2->getN();
    $k = 0;
    $ra = array_fill(0,$lv1*$lv2-1,0.0);

    for ($i = 0; $i < $lv1; $i++) {
        for ($j = 0; $j < $lv2; $j++) {
            $ra[$k] = $v1[$i]*$v2[$j];
            $k += 1;
        }
    }
    $r = new Vector($ra);
    return $r;
}

$cexp = function ($phi) { return new Complex(cos(phi),sin(phi)); };

$k1 = new Vector([0.0, 1.0]);
$k0 = new Vector([1.0, 0.0]);

$k00 = tensorProduct($k0,$k0);

$kk0 = new Matrix(
    [
        [1.0],
        [0.0]
    ]
);

$I = new Complex(0.0,1.0);
$MI = new Complex(0.0,-1.0);

$ck0 = new Vector([$I,$MI]);

$invsqr2 = 1.0/sqrt(2.0);

$H = new Matrix(
    [
        [1.0,  1.0],
        [1.0, -1.0]
    ]
);

$Y = new Matrix(
    [
        [0.0,  $I],
        [$MI, 0.0]
    ]
);

$R = function ($phi) { return new Matrix([[1.0,0.0],[0.0,$cexp($phi)]]);};

$H = $H->scalarMultiply($invsqr2);
$HH = $H->kroneckerProduct($H);

$kk00 = $kk0->kroneckerProduct($kk0);
$kkk00 = $kk00->asVectors();
$kkkk00 = new Vector($kkk00);

$prod  = $H->vectorMultiply($k0);
#$prod1 = $Y->vectorMultiply($k0);

function printq(Vector $v, int $width)
{
    $len = $v->getN();

    $noexpand = false;

    if ($width == 0) {
        $noexpand = true;
    }

    if ($noexpand) {
        for ($i = 0; $i < $len; $i++) {
            if ($v->get($i) <> 0) {
                echo $v->get($i) . "|" . decbin($i) . "> ";
            }
        }
    }

    echo "\n";
}

printq($prod,0);

#printq($prod1,0);

print $ck0;

$R(3.14);

