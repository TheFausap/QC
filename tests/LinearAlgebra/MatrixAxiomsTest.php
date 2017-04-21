<?php
namespace MathPHP\LinearAlgebra;

/**
 * Tests of Matrix axioms
 * These tests don't test specific functions,
 * but rather matrix axioms which in term make use of multiple functions.
 * If all the Matrix math is implemented properly, these tests should
 * all work out according to the axioms.
 *
 * Axioms tested:
 *  - Addition
 *    - r(A + B) = rA + rB
 *  - Multiplication
 *    - (AB)C = A(BC)
 *    - A(B + C) = AB + BC
 *    - r(AB) = (rA)B = A(rB)
 *  - Identity
 *    - AI = A = IA
 *  - Inverse
 *    - AA⁻¹ = I = A⁻¹A
 *    - (A⁻¹)⁻¹ = A
 *    - (AB)⁻¹ = B⁻¹A⁻¹
 *    - A is invertible, Aᵀ is inveritble
 *    - A is invertible, AAᵀ is inveritble
 *    - A is invertible, AᵀA is inveritble
 *  - Transpose
 *    - (Aᵀ)ᵀ = A
 *    - (A⁻¹)ᵀ = (Aᵀ)⁻¹
 *    - (rA)ᵀ = rAᵀ
 *    - (AB)ᵀ = BᵀAᵀ
 *    - (A + B)ᵀ = Aᵀ + Bᵀ
 *  - Trace
 *    - tr(A) = tr(Aᵀ)
 *    - tr(AB) = tr(BA)
 *  - Determinant
 *    - det(A) = det(Aᵀ)
 *    - det(AB) = det(A)det(B)
 *  - LU Decomposition (PA = LU)
 *    - PA = LU
 *    - A = P⁻¹LU
 *    - PPᵀ = I = PᵀP
 *    - (PA)⁻¹ = (LU)⁻¹ = U⁻¹L⁻¹
 *    - P⁻¹ = Pᵀ
 *  - System of linear equations (Ax = b)
 *    - Ax - b = 0
 *    - x = A⁻¹b
 *  - Symmetric matrix
 *    - A = Aᵀ
 *    - A⁻¹Aᵀ = I
 *    - A + B is symmetric
 *    - A - B is symmetric
 *    - kA is symmetric
 *    - AAᵀ is symmetric
 *    - AᵀA is symmetric
 *    - A is invertible symmetric, A⁻¹ is symmetric
 *  - Kronecker product
 *    - A ⊗ (B + C) = A ⊗ B + A ⊗ C
 *    - (A + B) ⊗ C = A ⊗ C + B ⊗ C
 *    - (A ⊗ B) ⊗ C = A ⊗ (B ⊗ C)
 *    - (kA) ⊗ B = A ⊗ (kB) = k(A ⊗ B)
 *    - (A ⊗ B)⁻¹ = A⁻¹ ⊗ B⁻¹
 *    - (A ⊗ B)ᵀ = Aᵀ ⊗ Bᵀ
 *    - det(A ⊗ B) = det(A)ᵐ det(B)ⁿ
 *  - Kronecker sum
 *    - A⊕B = A⊗Ib + I⊗aB
 *  - Covariance matrix
 *    - S = Sᵀ
 *  - Positive definiteness
 *    - A is PD ⇔ -A is ND
 *    - A is PSD ⇔ -A is NSD
 *    - A is PD ⇒ A is PSD
 *    - A is ND ⇒ A is NSD
 *    - A is PD ⇔ A⁻¹ is PD
 *    - A is ND ⇔ A⁻¹ is ND
 *    - A is PD and r > 0 ⇒ rA is PD
 *    - A and B are PD ⇒ A + B is PD
 *    - A and B are PD ⇒ ABA is PD
 *    - A and B are PD ⇒ BAB is PD
 *  - Triangular
 *    - Zero matrix is lower triangular
 *    - Zero matrix is upper triangular
 *    - Lᵀ is upper triangular
 *    - Uᵀ is lower triangular
 *    - LL is lower triangular
 *    - UU is upper triangular
 *    - L + L is lower triangular
 *    - U + U is upper triangular
 *    - L⁻¹ is lower triangular (If L is invertible)
 *    - U⁻¹ is upper triangular (If U is invertible)
 *    - kL is lower triangular
 *    - kU is upper triangular
 *    - L is invertible iif diagonal is all non zero
 *    - U is invertible iif diagonal is all non zero
 *  - Diagonal
 *    - Zero matrix is diagonal
 *    - Dᵀ is diagonal
 *    - DD is diagonal
 *    - D + D is diagonal
 *    - D⁻¹ is diagonal (If D is invertible)
 *    - kD is lower triangular
 *    - D is invertible iif diagonal is all non zero
 *  - Reduced row echelon form
 *    - RREF is upper triangular
 */
class MatrixAxiomsTest extends \PHPUnit_Framework_TestCase
{
    /**
     * Axiom: r(A + B) = rA + rB
     * Order of scalar multiplication does not matter.
     *
     * @dataProvider dataProviderForScalarMultiplicationOrderAddition
     */
    public function testScalarMultiplicationOrderAddition(array $A, array $B, int $r)
    {
        $A = MatrixFactory::create($A);
        $B = MatrixFactory::create($B);

        // r(A + B)
        $A＋B = $A->add($B);
        $r⟮A＋B⟯ = $A＋B->scalarMultiply($r);

        // rA + rB
        $rA     = $A->scalarMultiply($r);
        $rB     = $B->scalarMultiply($r);
        $rA＋rB = $rA->add($rB);

        $this->assertEquals($r⟮A＋B⟯->getMatrix(), $rA＋rB->getMatrix());
    }

    public function dataProviderForScalarMultiplicationOrderAddition()
    {
        return [
            [
                [
                    [1, 5],
                    [4, 3],
                ],
                [
                    [5, 6],
                    [2, 1],
                ], 5
            ],
            [
                [
                    [3, 8, 5],
                    [3, 6, 1],
                    [9, 5, 8],
                ],
                [
                    [5, 3, 8],
                    [6, 4, 5],
                    [1, 8, 9],
                ], 4
            ],
            [
                [
                    [-4, -2, 9],
                    [3, 14, -6],
                    [3, 9, 9],
                ],
                [
                    [8, 7, 8],
                    [-5, 4, 1],
                    [3, 5, 1],
                ], 7
            ],
            [
                [
                    [4, 7, 7, 8],
                    [3, 6, 4, 1],
                    [-3, 6, 8, -3],
                    [3, 2, 1, -54],
                ],
                [
                    [3, 2, 6, 7],
                    [4, 3, -6, 2],
                    [12, 14, 14, -6],
                    [4, 6, 4, -42],
                ], -8
            ],
        ];
    }

    /**
     * Axiom: (AB)C = A(BC)
     * Matrix multiplication is associative
     *
     * @dataProvider dataProviderForMultiplicationIsAssociative
     */
    public function testMultiplicationIsAssociative(array $A, array $B, array $C)
    {
        $A = MatrixFactory::create($A);
        $B = MatrixFactory::create($B);
        $C = MatrixFactory::create($C);

        $⟮AB⟯  = $A->multiply($B);
        $⟮AB⟯C = $⟮AB⟯->multiply($C);

        $⟮BC⟯  = $B->multiply($C);
        $A⟮BC⟯ = $A->multiply($⟮BC⟯);

        $this->assertEquals($⟮AB⟯C->getMatrix(), $A⟮BC⟯->getMatrix());
    }

    public function dataProviderForMultiplicationIsAssociative()
    {
        return [
            [
                [
                    [1, 5, 3],
                    [3, 6, 3],
                    [6, 7, 8],
                ],
                [
                    [6, 9, 9],
                    [3, 5, 1],
                    [3, 5, 12],
                ],
                [
                    [7, 4, 6],
                    [2, 3, 1],
                    [10, 12, 5],
                ],
            ],
            [
                [
                    [12, 21, 6],
                    [-3, 11, -6],
                    [3, 6, -3],
                ],
                [
                    [3, 7, 8],
                    [4, 4, 2],
                    [6, -4, 1],
                ],
                [
                    [1, -1, -5],
                    [6, 5, 19],
                    [3, 6, -2],
                ],
            ],
        ];
    }

    /**
     * Axiom: A(B + C) = AB + AC
     * Matrix multiplication is distributive
     *
     * @dataProvider dataProviderForMultiplicationIsDistributive
     */
    public function testMultiplicationIsDistributive(array $A, array $B, array $C)
    {
        $A = MatrixFactory::create($A);
        $B = MatrixFactory::create($B);
        $C = MatrixFactory::create($C);

        // A(B + C)
        $⟮B＋C⟯  = $B->add($C);
        $A⟮B＋C⟯ = $A->multiply($⟮B＋C⟯);

        // AB + AC
        $AB     = $A->multiply($B);
        $AC     = $A->multiply($C);
        $AB＋AC = $AB->add($AC);

        $this->assertEquals($A⟮B＋C⟯->getMatrix(), $AB＋AC->getMatrix());
    }

    public function dataProviderForMultiplicationIsDistributive()
    {
        return [
            [
                [
                    [1, 2],
                    [0, -1],
                ],
                [
                    [0, -1],
                    [1, 1],
                ],
                [
                    [-2, 0],
                    [0, 1],
                ],
            ],
            [
                [
                    [1, 5, 3],
                    [3, 6, 3],
                    [6, 7, 8],
                ],
                [
                    [6, 9, 9],
                    [3, 5, 1],
                    [3, 5, 12],
                ],
                [
                    [7, 4, 6],
                    [2, 3, 1],
                    [10, 12, 5],
                ],
            ],
            [
                [
                    [12, 21, 6],
                    [-3, 11, -6],
                    [3, 6, -3],
                ],
                [
                    [3, 7, 8],
                    [4, 4, 2],
                    [6, -4, 1],
                ],
                [
                    [1, -1, -5],
                    [6, 5, 19],
                    [3, 6, -2],
                ],
            ],
        ];
    }

    /**
     * Axiom: r(AB) = (rA)B = A(rB)
     * Order of scalar multiplication does not matter.
     *
     * @dataProvider dataProviderForScalarMultiplicationOrder
     */
    public function testScalarMultiplcationOrder(array $A, array $B, int $r)
    {
        $A = MatrixFactory::create($A);
        $B = MatrixFactory::create($B);

        // r(AB)
        $AB = $A->multiply($B);
        $r⟮AB⟯ = $AB->scalarMultiply($r);

        // (rA)B
        $rA = $A->scalarMultiply($r);
        $⟮rA⟯B = $rA->multiply($B);

        // A(rB)
        $rB = $B->scalarMultiply($r);
        $A⟮rB⟯ = $A->multiply($rB);

        $this->assertEquals($r⟮AB⟯->getMatrix(), $⟮rA⟯B->getMatrix());
        $this->assertEquals($⟮rA⟯B->getMatrix(), $A⟮rB⟯->getMatrix());
        $this->assertEquals($r⟮AB⟯->getMatrix(), $A⟮rB⟯->getMatrix());
    }

    public function dataProviderForScalarMultiplicationOrder()
    {
        return [
            [
                [
                    [1, 5],
                    [4, 3],
                ],
                [
                    [5, 6],
                    [2, 1],
                ], 5
            ],
            [
                [
                    [3, 8, 5],
                    [3, 6, 1],
                    [9, 5, 8],
                ],
                [
                    [5, 3, 8],
                    [6, 4, 5],
                    [1, 8, 9],
                ], 4
            ],
            [
                [
                    [-4, -2, 9],
                    [3, 14, -6],
                    [3, 9, 9],
                ],
                [
                    [8, 7, 8],
                    [-5, 4, 1],
                    [3, 5, 1],
                ], 7
            ],
            [
                [
                    [4, 7, 7, 8],
                    [3, 6, 4, 1],
                    [-3, 6, 8, -3],
                    [3, 2, 1, -54],
                ],
                [
                    [3, 2, 6, 7],
                    [4, 3, -6, 2],
                    [12, 14, 14, -6],
                    [4, 6, 4, -42],
                ], -8
            ],
        ];
    }

    /**
     * Axiom: AI = A = IA
     * Matrix multiplied with the identity matrix is the original matrix.
     *
     * @dataProvider dataProviderForMatrixTimesIdentityIsOriginalMatrix
     */
    public function testMatrixTimesIdentityIsOriginalMatrix(array $A)
    {
        $A  = MatrixFactory::create($A);
        $I  = MatrixFactory::identity($A->getN());
        $AI = $A->multiply($I);
        $IA = $I->multiply($A);

        $this->assertEquals($A->getMatrix(), $AI->getMatrix());
        $this->assertEquals($A->getMatrix(), $IA->getMatrix());
    }

    public function dataProviderForMatrixTimesIdentityIsOriginalMatrix()
    {
        return [
            [
                [
                    [1, 5],
                    [4, 3],
                ],
            ],
            [
                [
                    [5, 6],
                    [2, 1],
                ],
            ],
            [
                [
                    [3, 8, 5],
                    [3, 6, 1],
                    [9, 5, 8],
                ],
            ],
            [
                [
                    [5, 3, 8],
                    [6, 4, 5],
                    [1, 8, 9],
                ],
            ],
            [
                [
                    [-4, -2, 9],
                    [3, 14, -6],
                    [3, 9, 9],
                ],
            ],
            [
                [
                    [8, 7, 8],
                    [-5, 4, 1],
                    [3, 5, 1],
                ],
            ],
            [
                [
                    [4, 7, 7, 8],
                    [3, 6, 4, 1],
                    [-3, 6, 8, -3],
                    [3, 2, 1, -54],
                ],
            ],
            [
                [
                    [3, 2, 6, 7],
                    [4, 3, -6, 2],
                    [12, 14, 14, -6],
                    [4, 6, 4, -42],
                ],
            ],
        ];
    }

    /**
     * Axiom: AA⁻¹ = I = A⁻¹A
     * Matrix multiplied with its inverse is the identity matrix.
     *
     * @dataProvider dataProviderForInverse
     */
    public function testMatrixTimesInverseIsIdentity(array $A, array $A⁻¹)
    {
        $A    = MatrixFactory::create($A);
        $A⁻¹  = $A->inverse();
        $AA⁻¹ = $A->multiply($A⁻¹);
        $A⁻¹A = $A⁻¹->multiply($A);
        $I    = MatrixFactory::identity($A->getN());

        $this->assertEquals($I->getMatrix(), $AA⁻¹->getMatrix());
        $this->assertEquals($I->getMatrix(), $A⁻¹A->getMatrix());
        $this->assertEquals($AA⁻¹->getMatrix(), $A⁻¹A->getMatrix());
    }

    /**
     * Axiom: (A⁻¹)⁻¹ = A
     * Inverse of inverse is the original matrix.
     *
     * @dataProvider dataProviderForOneSquareMatrix
     */
    public function testInverseOfInverseIsOriginalMatrix(array $A)
    {
        $A       = MatrixFactory::create($A);
        $⟮A⁻¹⟯⁻¹  = $A->inverse()->inverse();

        $this->assertEquals($A->getMatrix(), $⟮A⁻¹⟯⁻¹->getMatrix());
    }

    public function dataProviderForInverse()
    {
        return [
            [
                [
                    [4, 7],
                    [2, 6],
                ],
                [
                    [0.6, -0.7],
                    [-0.2, 0.4],
                ],
            ],
            [
                [
                    [4, 3],
                    [3, 2],
                ],
                [
                    [-2, 3],
                    [3, -4],
                ],
            ],
            [
                [
                    [1, 2],
                    [3, 4],
                ],
                [
                    [-2, 1],
                    [3/2, -1/2],
                ],
            ],
            [
                [
                    [3, 3.5],
                    [3.2, 3.6],
                ],
                [
                    [-9, 8.75],
                    [8, -7.5],
                ],
            ],
            [
                [
                    [1, 2, 3],
                    [0, 4, 5],
                    [1, 0, 6],
                ],
                [
                    [12/11, -6/11, -1/11],
                    [5/22, 3/22, -5/22],
                    [-2/11, 1/11, 2/11],
                ],
            ],
            [
                [
                    [7, 2, 1],
                    [0, 3, -1],
                    [-3, 4, -2],
                ],
                [
                    [-2, 8, -5],
                    [3, -11, 7],
                    [9, -34, 21],
                ],
            ],
            [
                [
                    [3, 6, 6, 8],
                    [4, 5, 3, 2],
                    [2, 2, 2, 3],
                    [6, 8, 4, 2],
                ],
                [
                    [-0.333, 0.667, 0.667, -0.333],
                    [0.167, -2.333, 0.167, 1.417],
                    [0.167, 4.667, -1.833, -2.583],
                    [0.000, -2.000, 1.000, 1.000],
                ],
            ],
            [
                [
                    [4, 23, 6, 4, 7],
                    [3, 64, 23, 52, 2],
                    [65, 45, 3, 23, 1],
                    [2, 3, 4, 3, 9],
                    [53, 99, 54, 32, 105],
                ],
                [
                    [-0.142, 0.006, 0.003, -0.338, 0.038],
                    [0.172, -0.012, 0.010, 0.275, -0.035],
                    [-0.856, 0.082, -0.089, -2.344, 0.257],
                    [0.164, -0.001, 0.026, 0.683, -0.070],
                    [0.300, -0.033, 0.027, 0.909, -0.088],
                ],
            ],
        ];
    }

    /**
     * (AB)⁻¹ = B⁻¹A⁻¹
     * The inverse of a product is the reverse product of the inverses.
     *
     * @dataProvider dataProviderForInverseProductIsReverseProductOfInverses
     */
    public function testInverseProductIsReverseProductOfInverses(array $A, array $B)
    {
        $A = MatrixFactory::create($A);
        $B = MatrixFactory::create($B);

        $⟮AB⟯⁻¹ = $A->multiply($B)->inverse();

        $B⁻¹ = $B->inverse();
        $A⁻¹ = $A->inverse();
        $B⁻¹A⁻¹ = $B⁻¹->multiply($A⁻¹);

        $this->assertEquals($⟮AB⟯⁻¹->getMatrix(), $B⁻¹A⁻¹->getMatrix());
    }

    /**
     * @testCase Axiom: A is invertible, Aᵀ is inveritble
     * If A is an invertible matrix, then the tranpose is also inveritble
     * @dataProvider dataProviderForInverse
     * @param        array $A
     */
    public function testIfMatrixIsInvertibleThenTransposeIsInvertible(array $A)
    {
        $A  = MatrixFactory::create($A);
        $Aᵀ = $A->transpose();

        if ($A->isInvertible()) {
            $this->assertTrue($Aᵀ->isInvertible());
        } else {
            $this->assertFalse($Aᵀ->isInvertible());
        }
    }

    /**
     * @testCase Axiom: A is invertible, AAᵀ is inveritble
     * If A is an invertible matrix, then the product of A and tranpose of A is also inveritble
     * @dataProvider dataProviderForInverse
     * @param        array $A
     */
    public function testIfMatrixIsInvertibleThenProductOfMatrixAndTransposeIsInvertible(array $A)
    {
        $A   = MatrixFactory::create($A);
        $Aᵀ  = $A->transpose();
        $AAᵀ = $A->multiply($Aᵀ);

        if ($A->isInvertible()) {
            $this->assertTrue($AAᵀ->isInvertible());
        } else {
            $this->assertFalse($AAᵀ->isInvertible());
        }
    }

    /**
     * @testCase Axiom: A is invertible, AᵀA is inveritble
     * If A is an invertible matrix, then the product of transpose and A is also inveritble
     * @dataProvider dataProviderForInverse
     * @param        array $A
     */
    public function testIfMatrixIsInvertibleThenProductOfTransposeAndMatrixIsInvertible(array $A)
    {
        $A   = MatrixFactory::create($A);
        $Aᵀ  = $A->transpose();
        $AᵀA = $Aᵀ->multiply($A);

        if ($A->isInvertible()) {
            $this->assertTrue($AᵀA->isInvertible());
        } else {
            $this->assertFalse($AᵀA->isInvertible());
        }
    }

    public function dataProviderForInverseProductIsReverseProductOfInverses()
    {
        return [
            [
                [
                    [1, 5],
                    [4, 3],
                ],
                [
                    [5, 6],
                    [2, 1],
                ],
            ],
            [
                [
                    [3, 8, 5],
                    [3, 6, 1],
                    [9, 5, 8],
                ],
                [
                    [5, 3, 8],
                    [6, 4, 5],
                    [1, 8, 9],
                ],
            ],
            [
                [
                    [-4, -2, 9],
                    [3, 14, -6],
                    [3, 9, 9],
                ],
                [
                    [8, 7, 8],
                    [-5, 4, 1],
                    [3, 5, 1],
                ],
            ],
            [
                [
                    [4, 7, 7, 8],
                    [3, 6, 4, 1],
                    [-3, 6, 8, -3],
                    [3, 2, 1, -54],
                ],
                [
                    [3, 2, 6, 7],
                    [4, 3, -6, 2],
                    [12, 14, 14, -6],
                    [4, 6, 4, -42],
                ],
            ],
        ];
    }

    /**
     * (Aᵀ)ᵀ = A
     * The transpose of the transpose is the original matrix.
     *
     * @dataProvider dataProviderForOneSquareMatrix
     */
    public function testTransposeOfTransposeIsOriginalMatrix(array $A)
    {
        $A     = MatrixFactory::create($A);
        $⟮A⁻ᵀ⟯ᵀ = $A->transpose()->transpose();

        $this->assertEquals($⟮A⁻ᵀ⟯ᵀ->getMatrix(), $A->getMatrix());
    }

    /**
     * (A⁻¹)ᵀ = (Aᵀ)⁻¹
     * The transpose of the inverse is the inverse of the transpose.
     *
     * @dataProvider dataProviderForOneSquareMatrix
     */
    public function testTransposeOfInverseIsInverseOfTranspose(array $A)
    {
        $A     = MatrixFactory::create($A);
        $⟮A⁻¹⟯ᵀ = $A->inverse()->transpose();
        $⟮Aᵀ⟯⁻¹ = $A->transpose()->inverse();

        $this->assertEquals($⟮A⁻¹⟯ᵀ->getMatrix(), $⟮Aᵀ⟯⁻¹->getMatrix());
    }

    /**
     * (rA)ᵀ = rAᵀ
     * Scalar multiplication order does not matter for transpose
     *
     * @dataProvider dataProviderForOneSquareMatrix
     */
    public function testScalarMultiplicationOfTransposeOrder(array $A)
    {
        $A     = MatrixFactory::create($A);
        $r     = 4;

        $⟮rA⟯ᵀ = $A->scalarMultiply($r)->transpose();
        $rAᵀ  = $A->transpose()->scalarMultiply($r);

        $this->assertEquals($⟮rA⟯ᵀ->getMatrix(), $rAᵀ->getMatrix());
    }

    public function dataProviderForOneSquareMatrix()
    {
        return [
            [
                [
                    [1, 5],
                    [4, 3],
                ],
            ],
            [
                [
                    [5, 6],
                    [2, 1],
                ],
            ],
            [
                [
                    [3, 8, 5],
                    [3, 6, 1],
                    [9, 5, 8],
                ],
            ],
            [
                [
                    [5, 3, 8],
                    [6, 4, 5],
                    [1, 8, 9],
                ],
            ],
            [
                [
                    [-4, -2, 9],
                    [3, 14, -6],
                    [3, 9, 9],
                ],
            ],
            [
                [
                    [8, 7, 8],
                    [-5, 4, 1],
                    [3, 5, 1],
                ],
            ],
            [
                [
                    [4, 7, 7, 8],
                    [3, 6, 4, 1],
                    [-3, 6, 8, -3],
                    [3, 2, 1, -54],
                ],
            ],
            [
                [
                    [3, 2, 6, 7],
                    [4, 3, -6, 2],
                    [12, 14, 14, -6],
                    [4, 6, 4, -42],
                ],
            ],
        ];
    }

    /**
     * (AB)ᵀ = BᵀAᵀ
     * Transpose of a product of matrices equals the product of their transposes in reverse order.
     *
     * @dataProvider dataProviderForTwoSquareMatrices
     */
    public function testTransposeProductIsProductOfTranposesInReverseOrder(array $A, array $B)
    {
        $A = MatrixFactory::create($A);
        $B = MatrixFactory::create($B);

        // (AB)ᵀ
        $⟮AB⟯ᵀ = $A->multiply($B)->transpose();

        // BᵀAᵀ
        $Bᵀ   = $B->transpose();
        $Aᵀ   = $A->transpose();
        $BᵀAᵀ = $Bᵀ->multiply($Aᵀ);

        $this->assertEquals($⟮AB⟯ᵀ->getMatrix(), $BᵀAᵀ->getMatrix());
    }

    /**
     * (A + B)ᵀ = Aᵀ + Bᵀ
     * Transpose of sum is the same as sum of transposes
     *
     * @dataProvider dataProviderForTwoSquareMatrices
     */
    public function testTransposeSumIsSameAsSumOfTransposes(array $A, array $B)
    {
        $A = MatrixFactory::create($A);
        $B = MatrixFactory::create($B);

        // (A + B)ᵀ
        $⟮A＋B⟯ᵀ = $A->add($B)->transpose();

        // Aᵀ + Bᵀ
        $Aᵀ     = $A->transpose();
        $Bᵀ     = $B->transpose();
        $Aᵀ＋Bᵀ = $Aᵀ->add($Bᵀ);

        $this->assertEquals($⟮A＋B⟯ᵀ->getMatrix(), $Aᵀ＋Bᵀ->getMatrix());
    }

    public function dataProviderForTwoSquareMatrices()
    {
        return [
            [
                [
                    [1, 5],
                    [4, 3],
                ],
                [
                    [5, 6],
                    [2, 1],
                ],
            ],
            [
                [
                    [3, 8, 5],
                    [3, 6, 1],
                    [9, 5, 8],
                ],
                [
                    [5, 3, 8],
                    [6, 4, 5],
                    [1, 8, 9],
                ],
            ],
            [
                [
                    [-4, -2, 9],
                    [3, 14, -6],
                    [3, 9, 9],
                ],
                [
                    [8, 7, 8],
                    [-5, 4, 1],
                    [3, 5, 1],
                ],
            ],
            [
                [
                    [4, 7, 7, 8],
                    [3, 6, 4, 1],
                    [-3, 6, 8, -3],
                    [3, 2, 1, -54],
                ],
                [
                    [3, 2, 6, 7],
                    [4, 3, -6, 2],
                    [12, 14, 14, -6],
                    [4, 6, 4, -42],
                ],
            ],
        ];
    }

    /**
     * tr(A) = tr(Aᵀ)
     * Trace is the same as the trace of the transpose
     *
     * @dataProvider dataProviderForOneSquareMatrix
     */
    public function testTraceIsSameAsTraceOfTranspose(array $A)
    {
        $A  = MatrixFactory::create($A);
        $Aᵀ = $A->transpose();

        $tr⟮A⟯  = $A->trace();
        $tr⟮Aᵀ⟯ = $Aᵀ->trace();

        $this->assertEquals($tr⟮A⟯, $tr⟮Aᵀ⟯);
    }

    /**
     * tr(AB) = tr(BA)
     * Trace of product does not matter the order they were multiplied
     *
     * @dataProvider dataProviderForTwoSquareMatrices
     */
    public function testTraceOfProductIsSameRegardlessOfOrderMultiplied(array $A, array $B)
    {
        $A  = MatrixFactory::create($A);
        $B  = MatrixFactory::create($B);

        $tr⟮AB⟯ = $A->multiply($B)->trace();
        $tr⟮BA⟯ = $B->multiply($A)->trace();

        $this->assertEquals($tr⟮AB⟯, $tr⟮BA⟯);
    }

    /**
     * det(A) = det(Aᵀ)
     * Determinant of matrix is the same as determinant of transpose.
     *
     * @dataProvider dataProviderForOneSquareMatrix
     */
    public function testDeterminantSameAsDeterminantOfTranspose(array $A)
    {
        $A = MatrixFactory::create($A);

        $det⟮A⟯  = $A->det();
        $det⟮Aᵀ⟯ = $A->transpose()->det();

        $this->assertEquals($det⟮A⟯, $det⟮Aᵀ⟯);
    }

    /**
     * det(AB) = det(A)det(B)
     * Determinant of product of matrices is the same as the product of determinants.
     *
     * @dataProvider dataProviderForTwoSquareMatrices
     */
    public function testDeterminantProductSameAsProductOfDeterminants(array $A, array $B)
    {
        $A = MatrixFactory::create($A);
        $B = MatrixFactory::create($B);

        // det(AB)
        $det⟮AB⟯  = $A->multiply($B)->det();

        // det(A)det(B)
        $det⟮A⟯ = $A->det();
        $det⟮B⟯ = $B->det();
        $det⟮A⟯det⟮B⟯ = $det⟮A⟯ * $det⟮B⟯;

        $this->assertEquals($det⟮AB⟯, $det⟮A⟯det⟮B⟯, '', 0.000001);
    }

    /**
     * PA = LU
     * Basic LU decomposition property that permutation matrix times the matrix is the product of the lower and upper decomposition matrices.
     * @dataProvider dataProviderForOneSquareMatrix
     */
    public function testLUDecompositionPAEqualsLU(array $A)
    {
        $A = MatrixFactory::create($A);

        $LUP = $A->LUDecomposition();

        $L   = $LUP['L'];
        $U   = $LUP['U'];
        $P   = $LUP['P'];

        $PA = $P->multiply($A);
        $LU = $L->multiply($U);
        $UL = $U->multiply($L);

        $this->assertEquals($PA->getMatrix(), $LU->getMatrix());
    }

    /**
     * A = P⁻¹LU
     * @dataProvider dataProviderForOneSquareMatrix
     */
    public function testLUDecompositionAEqualsPInverseLU(array $A)
    {
        $A = MatrixFactory::create($A);

        $LUP = $A->LUDecomposition();

        $L   = $LUP['L'];
        $U   = $LUP['U'];
        $P   = $LUP['P'];

        $P⁻¹LU = $P->inverse()->multiply($L)->multiply($U);

        $this->assertEquals($A->getMatrix(), $P⁻¹LU->getMatrix());
    }

    /**
     * PPᵀ = I = PᵀP
     * Permutation matrix of the LU decomposition times the transpose of the permutation matrix is the identity matrix.
     * @dataProvider dataProviderForOneSquareMatrix
     */
    public function testLUDecompositionPPTransposeEqualsIdentity(array $A)
    {
        $A = MatrixFactory::create($A);

        $LUP = $A->LUDecomposition();

        $P  = $LUP['P'];
        $Pᵀ = $P->transpose();

        $PPᵀ = $P->multiply($Pᵀ);
        $PᵀP = $Pᵀ->multiply($P);

        $I = MatrixFactory::identity($A->getM());

        $this->assertEquals($PPᵀ->getMatrix(), $I->getMatrix());
        $this->assertEquals($I->getMatrix(), $PᵀP->getMatrix());
        $this->assertEquals($PPᵀ->getMatrix(), $PᵀP->getMatrix());
    }

    /**
     * (PA)⁻¹ = (LU)⁻¹ = U⁻¹L⁻¹
     * Inverse of the LU decomposition equation is the inverse of the other side.
     * @dataProvider dataProviderForOneSquareMatrix
     */
    public function testInverseWithLUDecompositionInverse(array $A)
    {
        $A = MatrixFactory::create($A);

        $LUP = $A->LUDecomposition();
        $L   = $LUP['L'];
        $U   = $LUP['U'];
        $P   = $LUP['P'];

        $⟮PA⟯⁻¹    = $P->multiply($A)->inverse();
        $⟮LU⟯⁻¹ = $L->multiply($U)->inverse();

        $U⁻¹      = $U->inverse();
        $L⁻¹      = $L->inverse();
        $U⁻¹L⁻¹   = $U⁻¹->multiply($L⁻¹);

        $this->assertEquals($⟮PA⟯⁻¹->getMatrix(), $⟮LU⟯⁻¹->getMatrix());
        $this->assertEquals($⟮LU⟯⁻¹->getMatrix(), $U⁻¹L⁻¹->getMatrix());
        $this->assertEquals($⟮PA⟯⁻¹->getMatrix(), $U⁻¹L⁻¹->getMatrix());
    }

    /**
     * P⁻¹ = Pᵀ
     * Inverse of the permutation matrix equals the transpose of the permutation matrix
     * @dataProvider dataProviderForOneSquareMatrix
     */
    public function testPInverseEqualsPTranspose(array $A)
    {
        $A = MatrixFactory::create($A);

        $LUP = $A->LUDecomposition();
        $P   = $LUP['P'];
        $P⁻¹ = $P->inverse();
        $Pᵀ  = $P->transpose();

        $this->assertEquals($P⁻¹, $Pᵀ);
    }

    /**
     * Axiom: Ax - b = 0
     * Matrix multiplied with unknown x vector subtract solution b is 0
     *
     * @dataProvider dataProviderForSolve
     */
    public function testSolveEquationForZero(array $A, array $b, array $x, array $zeros)
    {
        $A = MatrixFactory::create($A);
        $x = new Vector($x);
        $b = (new Vector($b))->asColumnMatrix();
        $z = (new Vector($zeros))->asColumnMatrix();

        // Ax - b
        $R = $A->multiply($x)->subtract($b);

        $this->assertEquals($z, $R, '', 0.01);
    }

    /**
     * Axiom: x = A⁻¹b
     * Matrix multiplied with unknown x vector subtract solution b is 0
     *
     * @dataProvider dataProviderForSolve
     */
    public function testSolveInverseBEqualsX(array $A, array $b, array $x, array $zeros)
    {
        $A   = MatrixFactory::create($A);
        $A⁻¹ = $A->inverse();
        $x = (new Vector($x))->asColumnMatrix();
        $b = new Vector($b);

        // A⁻¹b
        $A⁻¹b = $A⁻¹->multiply($b);

        $this->assertEquals($x, $A⁻¹b, '', 0.001);
    }

    public function dataProviderForSolve()
    {
        return [
            [
                [
                    [3, 4],
                    [2, -1],
                ],
                [5, 7],
                [3, -1],
                [0, 0],
            ],
            [
                [
                    [3, 1],
                    [2, -1],
                ],
                [5, 0],
                [1, 2],
                [0, 0],
            ],
            [
                [
                    [3, 4],
                    [5, 3],
                ],
                [-2, 4],
                [2, -2],
                [0, 0],
            ],
            [
                [
                    [1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1],
                ],
                [2, 3, -4],
                [2, 3, -4],
                [0, 0, 0],
            ],
            [
                [
                    [1, 1, -1],
                    [3, 1, 1],
                    [1, -1, 4],
                ],
                [1, 9, 8],
                [3, -1, 1],
                [0, 0, 0],
            ],
            [
                [
                    [2, 4, 1],
                    [4, -10, 2],
                    [1, 2, 4],
                ],
                [5, -8, 13],
                [-1, 1, 3],
                [0, 0, 0],
            ],
            [
                [
                    [1, 1, 1],
                    [0, 2, 5],
                    [2, 5, -1],
                ],
                [6, -4, 27],
                [5, 3, -2],
                [0, 0, 0],
            ],
            [
                [
                    [1, 2, 3],
                    [2, -1, 1],
                    [3, 0, -1],
                ],
                [9, 8, 3],
                [2, -1, 3],
                [0, 0, 0],
            ],
            [
                [
                    [2, 1, -3],
                    [4, -2, 1],
                    [3, 5, -2],
                ],
                [-4, 9, 5],
                [2, 1, 3],
                [0, 0, 0],
            ],
            [
                [
                    [4, 9, 0],
                    [8, 0, 6],
                    [0, 6, 6],
                ],
                [8, -1, -1],
                [1/2, 2/3, -5/6],
                [0, 0, 0],
            ],
            [
                [
                    [1, 1, 1],
                    [1, -2, 2],
                    [1, 2, -1],
                ],
                [0, 4, 2],
                [4, -2, -2],
                [0, 0, 0],
            ],
            [
                [
                    [3, 3, 4],
                    [3, 5, 9],
                    [5, 9, 17],
                ],
                [1, 2, 4],
                [1, -2, 1],
                [0, 0, 0],
            ],
            [
                [
                    [2, 1, 1],
                    [-1, 1, -1],
                    [1, 2, 3],
                ],
                [2, 3, -10],
                [3, 1, -5],
                [0, 0, 0],
            ],
            [
                [
                    [4, 2, -1, 3],
                    [3, -4, 2, 5],
                    [-2, 6, -5, -2],
                    [5, 1, 6, -3],
                ],
                [16.9, -14, 25, 9.4],
                [4.5, 1.6, -3.8, -2.7],
                [0, 0, 0, 0],
            ],
            [
                [
                    [4, 2, -1, 3],
                    [3, -4, 2, 5],
                    [-2, 6, -5, -2],
                    [5, 1, 6, -3],
                ],
                [-12, 34, 27, -19],
                [-101.485, 101.242, 115.727, 102.394],
                [0, 0, 0, 0],
            ],
            [
                [
                    [ 4,  1,  2,  -3],
                    [-3,  3, -1,   4],
                    [-1,  2,  5,   1],
                    [ 5,  4,  3,  -1],
                ],
                [-16, 20, -4, -10],
                [-1, 1, -2, 3],
                [0, 0, 0, 0],
            ],
            [
                [
                    [ 4,  1,  2,  -3,  5],
                    [-3,  3, -1,   4, -2],
                    [-1,  2,  5,   1,  3],
                    [ 5,  4,  3,  -1,  2],
                    [ 1, -2,  3,  -4,  5],
                ],
                [-16, 20, -4, -10,  3],
                [-15.354, 15.813, -1.770, -22.148, -6.660],
                [0, 0, 0, 0, 0],
            ],
            [
                [
                    [1, 1, -2, 1, 3, -1],
                    [2, -1, 1, 2, 1, -3],
                    [1, 3, -3, -1, 2, 1],
                    [5, 2, -1, -1, 2, 1],
                    [-3, -1, 2, 3, 1, 3],
                    [4, 3, 1, -6, -3, -2],
                ],
                [4, 20, -15, -3, 16, -27],
                [1, -2, 3, 4, 2, -1],
                [0, 0, 0, 0, 0, 0],
            ],
        ];
    }

    /**
     * Axiom: A = Aᵀ
     * Symmetric matrix is the same as its transpose
     * @dataProvider dataProviderForSymmetric
     */
    public function testSymmetricEqualsTranspose(array $A)
    {
        $A  = MatrixFactory::create($A);
        $Aᵀ = $A->transpose();

        $this->assertEquals($A, $Aᵀ);
        $this->assertEquals($A->getMatrix(), $Aᵀ->getMatrix());
    }

    /**
     * Axiom: A⁻¹Aᵀ = I
     * Symmetric matrix inverse times tranpose equals identity matrix
     * @dataProvider dataProviderForSymmetric
     */
    public function testSymmetricInverseTranposeEqualsIdentity(array $A)
    {
        $A   = MatrixFactory::create($A);
        $A⁻¹ = $A->inverse();
        $Aᵀ  = $A->transpose();

        $A⁻¹Aᵀ = $A⁻¹->multiply($Aᵀ);
        $I     = MatrixFactory::identity($A->getM());

        $this->assertEquals($I, $A⁻¹Aᵀ);
        $this->assertEquals($I->getMatrix(), $A⁻¹Aᵀ->getMatrix());
    }

    /**
     * @testCase Axiom: A + B is symmetric
     * If A and B are symmetric matrices with the sme size, then A + B is symmetric
     * @dataProvider dataProviderForSymmetric
     * @param array $A
     */
    public function testSymmetricMatricesSumIsSymmetric(array $M)
    {
        $A    = MatrixFactory::create($M);
        $B    = MatrixFactory::create($M);
        $A＋B = $A->add($B);

        $this->assertTrue($A->isSymmetric());
        $this->assertTrue($B->isSymmetric());
        $this->assertTrue($A＋B->isSymmetric());
    }

    /**
     * @testCase Axiom: A - B is symmetric
     * If A and B are symmetric matrices with the sme size, then A - B is symmetric
     * @dataProvider dataProviderForSymmetric
     * @param array $A
     */
    public function testSymmetricMatricesDifferenceIsSymmetric(array $M)
    {
        $A   = MatrixFactory::create($M);
        $B   = MatrixFactory::create($M);
        $A−B = $A->subtract($B);

        $this->assertTrue($A->isSymmetric());
        $this->assertTrue($B->isSymmetric());
        $this->assertTrue($A−B->isSymmetric());
    }

    /**
     * @testCase Axiom: kA is symmetric
     * If A is a symmetric matrix, kA is symmetric
     * @dataProvider dataProviderForSymmetric
     * @param array $A
     */
    public function testSymmetricMatricesTimesScalarIsSymmetric(array $M)
    {
        $A   = MatrixFactory::create($M);
        $this->assertTrue($A->isSymmetric());

        foreach (range(1, 10) as $k) {
            $kA = $A->scalarMultiply($k);
            $this->assertTrue($kA->isSymmetric());
        }
    }

    /**
     * @testCase Axiom: AAᵀ is symmetric
     * If A is a symmetric matrix, AAᵀ is symmetric
     * @dataProvider dataProviderForSymmetric
     * @param array $A
     */
    public function testSymmetricMatrixTimesTransposeIsSymmetric(array $M)
    {
        $A   = MatrixFactory::create($M);
        $Aᵀ  = $A->transpose();
        $AAᵀ = $A->multiply($Aᵀ);

        $this->assertTrue($A->isSymmetric());
        $this->assertTrue($AAᵀ->isSymmetric());
    }

    /**
     * @testCase Axiom: AᵀA is symmetric
     * If A is a symmetric matrix, AᵀA is symmetric
     * @dataProvider dataProviderForSymmetric
     * @param array $A
     */
    public function testTransposeTimesSymmetricMatrixIsSymmetric(array $M)
    {
        $A   = MatrixFactory::create($M);
        $Aᵀ  = $A->transpose();
        $AᵀA = $Aᵀ->multiply($A);

        $this->assertTrue($A->isSymmetric());
        $this->assertTrue($AᵀA->isSymmetric());
    }

    /**
     * @testCase Axiom: A is invertible symmetric, A⁻¹ is symmetric
     * If A is an invertible symmetric matrix, the inverse of A is also symmetric
     * @dataProvider dataProviderForSymmetric
     * @param array $A
     */
    public function testMatrixIsInvertibleSummetricThenInverseIsSymmetric(array $M)
    {
        $A   = MatrixFactory::create($M);

        if ($A->isInvertible() && $A->isSymmetric()) {
            $A⁻¹ = $A->inverse();
            $A⁻¹ = $A->map(
                function ($x) {
                    return round($x, 5); // Floating point adjustment
                }
            );
            $this->assertTrue($A⁻¹->isSymmetric());
        }
    }

    public function dataProviderForSymmetric()
    {
        return [
            [
                [
                    [1],
                ],
            ],
            [
                [
                    [1, 2],
                    [2, 1],
                ],
            ],
            [
                [
                    [4, 1],
                    [1, -2],
                ],
            ],
            [
                [
                    [4, -1],
                    [-1, 9],
                ],
            ],
            [
                [
                    [1, 2, 3],
                    [2, 6, 4],
                    [3, 4, 5],
                ],
            ],
            [
                [
                    [1, 7, 3],
                    [7, 4, -5],
                    [3, -5, 6],
                ],
            ],
            [
                [
                    [5, 6, 7],
                    [6, 3, 2],
                    [7, 2, 1],
                ],
            ],
            [
                [
                    [2, 7, 3],
                    [7, 9, 4],
                    [3, 4, 7],
                ],
            ],
            [
                [
                    [4, -1, -1, -1],
                    [-1, 4, -1, -1],
                    [-1, -1, 4, -1],
                    [-1, -1, -1, 4],
                ],
            ],
            [
                [
                    [1, 5, 6, 8],
                    [5, 2, 7, 9],
                    [6, 7, 3, 10],
                    [8, 9, 10, 4],
                ],
            ],
        ];
    }

    /**
     * Axiom: A ⊗ (B + C) = A ⊗ B + A ⊗ C
     * Kronecker product bilinearity
     * @dataProvider dataProviderForThreeMatrices
     */
    public function testKroneckerProductBilinearity1(array $A, array $B, array $C)
    {
        $A = MatrixFactory::create($A);
        $B = MatrixFactory::create($B);
        $C = MatrixFactory::create($C);

        $A⊗⟮B＋C⟯  = $A->kroneckerProduct($B->add($C));
        $A⊗B＋A⊗C = $A->kroneckerProduct($B)->add($A->kroneckerProduct($C));

        $this->assertEquals($A⊗⟮B＋C⟯->getMatrix(), $A⊗B＋A⊗C->getMatrix());
    }

    /**
     * Axiom: (A + B) ⊗ C = A ⊗ C + B ⊗ C
     * Kronecker product bilinearity
     * @dataProvider dataProviderForThreeMatrices
     */
    public function testKroneckerProductBilinearity2(array $A, array $B, array $C)
    {
        $A = MatrixFactory::create($A);
        $B = MatrixFactory::create($B);
        $C = MatrixFactory::create($C);

        $⟮A＋B⟯⊗C  = $A->add($B)->kroneckerProduct($C);
        $A⊗C＋B⊗C = $A->kroneckerProduct($C)->add($B->kroneckerProduct($C));

        $this->assertEquals($⟮A＋B⟯⊗C->getMatrix(), $A⊗C＋B⊗C->getMatrix());
    }

    /**
     * Axiom: (A ⊗ B) ⊗ C = A ⊗ (B ⊗ C)
     * Kronecker product associative
     * @dataProvider dataProviderForThreeMatrices
     */
    public function testKroneckerProductAssociativity(array $A, array $B, array $C)
    {
        $A = MatrixFactory::create($A);
        $B = MatrixFactory::create($B);
        $C = MatrixFactory::create($C);

        $⟮A⊗B⟯⊗C = $A->kroneckerProduct($B)->kroneckerProduct($C);
        $A⊗⟮B⊗C⟯ = $A->kroneckerProduct($B->kroneckerProduct($C));

        $this->assertEquals($⟮A⊗B⟯⊗C->getMatrix(), $A⊗⟮B⊗C⟯->getMatrix());
    }

    public function dataProviderForThreeMatrices()
    {
        return [
            [
                [
                    [1],
                ],
                [
                    [2],
                ],
                [
                    [3],
                ],
            ],
            [
                [
                    [1, 5, 3],
                    [3, 6, 3],
                    [6, 7, 8],
                ],
                [
                    [6, 9, 9],
                    [3, 5, 1],
                    [3, 5, 12],
                ],
                [
                    [7, 9, 6],
                    [1, 9, 1],
                    [10, 12, 4],
                ],
            ],
            [
                [
                    [12, 21, 6],
                    [-3, 11, -6],
                    [3, 6, -3],
                ],
                [
                    [3, 7, 8],
                    [4, 4, 2],
                    [6, -4, 1],
                ],
                [
                    [-1, -1, -5],
                    [8, 15, 15],
                    [8, 6, -12],
                ],
            ],
            [
                [
                    [1, 2],
                    [0, -1],
                ],
                [
                    [0, -1],
                    [1, 1],
                ],
                [
                    [2, 8],
                    [2, 1],
                ],
            ],
            [
                [
                    [1, 5, 3],
                    [3, 6, 3],
                    [6, 7, 8],
                ],
                [
                    [6, 9, 9],
                    [3, 5, 1],
                    [3, 5, 12],
                ],
                [
                    [6, 4, 9],
                    [12, 3, -1],
                    [10, 2, 15],
                ],
            ],
            [
                [
                    [12, 21, 6],
                    [-3, 11, -6],
                    [3, 6, -3],
                ],
                [
                    [3, 7, 8],
                    [4, 4, 2],
                    [6, -4, 1],
                ],
                [
                    [1, 1, 5],
                    [3, 4, 9],
                    [3, 16, -2],
                ],
            ],
            [
                [
                    [1, 2, 3, 4, 5],
                    [2, 3, 4, 5, 6],
                    [4, 5, 6, 7, 8],
                    [6, 5, 4, 5, 7],
                ],
                [
                    [1, 2, 5, 5, 6],
                    [2, 3, 5, 5, 6],
                    [5, 4, 5, 5, 6],
                    [3, 2, 5, 5, 6],
                ],
                [
                    [5, 5, 7, 8, 9],
                    [4, 4, 7, 8, 9],
                    [7, 6, 7, 6, 7],
                    [9, 9, 9, 0, 0],
                ]
            ],
        ];
    }

    /**
     * Axiom: (kA) ⊗ B = A ⊗ (kB) = k(A ⊗ B)
     * Kronecker product scalar multiplication
     * @dataProvider dataProviderForTwoSquareMatrices
     */
    public function testKroneckerProductScalarMultiplication(array $A, array $B)
    {
        $A = MatrixFactory::create($A);
        $B = MatrixFactory::create($B);
        $k = 5;

        $⟮kA⟯⊗B = $A->scalarMultiply($k)->kroneckerProduct($B);
        $A⊗⟮kB⟯ = $A->kroneckerProduct($B->scalarMultiply($k));
        $k⟮A⊗B⟯ = $A->kroneckerProduct($B)->scalarMultiply($k);

        $this->assertEquals($⟮kA⟯⊗B->getMatrix(), $A⊗⟮kB⟯->getMatrix());
        $this->assertEquals($⟮kA⟯⊗B->getMatrix(), $k⟮A⊗B⟯->getMatrix());
        $this->assertEquals($k⟮A⊗B⟯->getMatrix(), $A⊗⟮kB⟯->getMatrix());
    }

    /**
     * Axiom: (A ⊗ B)⁻¹ = A⁻¹ ⊗ B⁻¹
     * Inverse of Kronecker product
     * @dataProvider dataProviderForTwoSquareMatrices
     */
    public function testKroneckerProductInverse(array $A, array $B)
    {
        $A = MatrixFactory::create($A);
        $B = MatrixFactory::create($B);

        $A⁻¹     = $A->inverse();
        $B⁻¹     = $B->inverse();
        $A⁻¹⊗B⁻¹ = $A⁻¹->kroneckerProduct($B⁻¹);
        $⟮A⊗B⟯⁻¹  = $A->kroneckerProduct($B)->inverse();

        $this->assertEquals($A⁻¹⊗B⁻¹->getMatrix(), $⟮A⊗B⟯⁻¹->getMatrix());
    }

    /**
     * Axiom: (A ⊗ B)ᵀ = Aᵀ ⊗ Bᵀ
     * Transpose of Kronecker product
     * @dataProvider dataProviderForTwoSquareMatrices
     */
    public function testKroneckerProductTranspose(array $A, array $B)
    {
        $A = MatrixFactory::create($A);
        $B = MatrixFactory::create($B);

        $Aᵀ    = $A->transpose();
        $Bᵀ    = $B->transpose();
        $Aᵀ⊗Bᵀ = $Aᵀ->kroneckerProduct($Bᵀ);
        $⟮A⊗B⟯ᵀ = $A->kroneckerProduct($B)->transpose();

        $this->assertEquals($Aᵀ⊗Bᵀ->getMatrix(), $⟮A⊗B⟯ᵀ->getMatrix());
    }

    /**
     * Axiom: det(A ⊗ B) = det(A)ᵐ det(B)ⁿ
     * Determinant of Kronecker product - where A is nxn matrix, and b is nxn matrix
     * @dataProvider dataProviderForKroneckerProductDeterminant
     */
    public function testKroneckerProductDeterminant(array $A, array $B)
    {
        $A = MatrixFactory::create($A);
        $B = MatrixFactory::create($B);

        $det⟮A⟯ᵐ  = ($A->det())**$B->getM();
        $det⟮B⟯ⁿ  = ($B->det())**$A->getN();
        $det⟮A⊗B⟯ = $A->kroneckerProduct($B)->det();

        $this->assertEquals($det⟮A⊗B⟯, $det⟮A⟯ᵐ  * $det⟮B⟯ⁿ, '', 0.0001);
    }

    public function dataProviderForKroneckerProductDeterminant()
    {
        return [
            [
                [
                    [5],
                ],
                [
                    [4],
                ],
            ],
            [
                [
                    [5, 6],
                    [2, 4],
                ],
                [
                    [4, 9],
                    [3, 1],
                ],
            ],
            [
                [
                    [5, 6],
                    [-2, 4],
                ],
                [
                    [4, -9],
                    [3, 1],
                ],
            ],
            [
                [
                    [1, 2, 3],
                    [2, 4, 6],
                    [7, 6, 5],
                ],
                [
                    [2, 3, 4],
                    [3, 1, 6],
                    [4, 3, 3],
                ],
            ],
            [
                [
                    [1, -2, 3],
                    [2, -4, 6],
                    [7, -6, 5],
                ],
                [
                    [2, 3, 4],
                    [3, 1, 6],
                    [4, 3, 3],
                ],
            ],
            [
                [
                    [1, 2, 3, 4],
                    [2, 4, 6, 8],
                    [7, 6, 5, 4],
                    [1, 3, 5, 7],
                ],
                [
                    [2, 3, 4, 5],
                    [3, 1, 6, 3],
                    [4, 3, 3, 4],
                    [3, 3, 4, 1],
                ],
            ],
            [
                [
                    [-1, 2, 3, 4],
                    [2, 4, 6, 8],
                    [7, 6, 5, 4],
                    [1, 3, 5, 7],
                ],
                [
                    [2, 3, 4, 5],
                    [3, 1, 6, 3],
                    [4, 3, 3, -4],
                    [3, 3, 4, 1],
                ],
            ],
            [
                [
                    [1, 2, 3, 4, 5],
                    [2, 4, 6, 8, 10],
                    [7, 6, 5, 4, 5],
                    [1, 3, 5, 7, 9],
                    [1, 2, 3, 4, 5],
                ],
                [
                    [2, 3, 4, 5, 2],
                    [3, 1, 6, 3, 2],
                    [4, 3, 3, 4, 2],
                    [3, 3, 4, 1, 2],
                    [1, 1, 2, 2, 1],
                ],
            ],
        ];
    }

    /**
     * Axiom: Kronecker sum A⊕B = A⊗Ib + I⊗aB
     * Kronecker sum is the matrix product of the Kronecker product of each matrix with the other matrix's identiry matrix.
     * @dataProvider dataProviderForTwoSquareMatrices
     */
    public function testKroneckerSum(array $A, array $B)
    {
        $A   = new SquareMatrix($A);
        $B   = new SquareMatrix($B);
        $A⊕B = $A->kroneckerSum($B);


        $In          = MatrixFactory::identity($A->getN());
        $Im          = MatrixFactory::identity($B->getM());
        $A⊗Im        = $A->kroneckerProduct($Im);
        $In⊗B        = $In->kroneckerProduct($B);
        $A⊗Im＋In⊗B = $A⊗Im->add($In⊗B);

        $this->assertEquals($A⊕B, $A⊗Im＋In⊗B);
    }

    /**
     * Axiom: Covariance matrix S = Sᵀ
     * Covariance matrix is symmetric so it is the same as its transpose
     * @dataProvider dataProviderForCovarianceSymmetric
     * @param        array $A
     */
    public function testCovarianceMatrixIsSymmetric(array $A)
    {
        $A  = MatrixFactory::create($A);
        $S  = $A->covarianceMatrix();
        $Sᵀ = $S->transpose();

        $this->assertEquals($S, $Sᵀ);
        $this->assertEquals($S->getMatrix(), $Sᵀ->getMatrix());
    }

    public function dataProviderForCovarianceSymmetric()
    {
        return [
            [
                [
                    [1, 4, 7, 8],
                    [2, 2, 8, 4],
                    [1, 13, 1, 5],
                ],
            ],
            [
                [
                    [19, 22, 6, 3, 2, 20],
                    [12, 6, 9, 15, 13, 5],
                ],
            ],
            [
                [
                    [4, 4.2, 3.9, 4.3, 4.1],
                    [2, 2.1, 2, 2.1, 2.2],
                    [.6, .59, .58, .62, .63]
                ],
            ],
            [
                [
                    [2.5, 0.5, 2.2, 1.9, 3.1, 2.3, 2, 1, 1.5, 1.1],
                    [2.4, 0.7, 2.9, 2.2, 3.0, 2.7, 1.6, 1.1, 1.6, 0.9],
                ],
            ],
        ];
    }

    /**
     * @testCase Axiom: Positive definiteness A is PD ⇔ -A is ND
     * If A is positive definite, then -A is negative definite.
     * @dataProvider dataProviderForPositiveDefiniteMatrix
     * @param        array $A
     */
    public function testPositiveDefiniteNegativeisNegativeDefinite(array $A)
    {
        $A  = MatrixFactory::create($A);
        $−A = $A->scalarMultiply(-1);

        $this->assertTrue($A->isPositiveDefinite());
        $this->assertTrue($−A->isNegativeDefinite());
    }

    public function dataProviderForPositiveDefiniteMatrix(): array
    {
        return [
            [
                [
                    [2, -1],
                    [-1, 2],
                ],
            ],
            [
                [
                    [1, -1],
                    [-1, 4],
                ],
            ],
            [
                [
                    [5, 2],
                    [2, 3],
                ],
            ],
            [
                [
                    [6, 4],
                    [4, 5],
                ],
            ],
            [
                [
                    [12, -12],
                    [-12, 96],
                ],
            ],
            [
                [
                    [2, -1, 0],
                    [-1, 2, -1],
                    [0, -1, 2],
                ],
            ],
            [
                [
                    [2, -1, 1],
                    [-1, 2, -1],
                    [1, -1, 2],
                ],
            ],
            [
                [
                    [1, 0, 0],
                    [0, 3, 0],
                    [0, 0, 2],
                ],
            ],
            [
                [
                    [3, -2, 0],
                    [-2, 2, 0],
                    [0, 0, 2],
                ],
            ],
            [
                [
                    [4, 1, -1],
                    [1, 2, 1],
                    [-1, 1, 2],
                ],
            ],
            [
                [
                    [9, -3, 3, 9],
                    [-3, 17, -1, -7],
                    [3, -1, 17, 15],
                    [9, -7, 15, 44],
                ],
            ],
            [
                [
                    [14, 4, 9],
                    [4, 14, -7],
                    [9, -7, 14],
                ],
            ],
            [
                [
                    [13, 0, -3],
                    [0, 9, 9],
                    [-3, 9, 10],
                ],
            ],
            [
                [
                    [14, -7, -13],
                    [-7, 6, 5],
                    [-13, 5, 14],
                ],
            ],
        ];
    }

    /**
     * @testCase Axiom: Positive semidefiniteness A is PSD ⇔ -A is NSD
     * If A is positive semidefinite, then -A is negative definite.
     * @dataProvider dataProviderForPositiveSemidefiniteMatrix
     * @param        array $A
     */
    public function testPositiveSemidefiniteNegativeisNegativeSemidefinite(array $A)
    {
        $A  = MatrixFactory::create($A);
        $−A = $A->scalarMultiply(-1);

        $this->assertTrue($A->isPositiveSemidefinite());
        $this->assertTrue($−A->isNegativeSemidefinite());
    }

    public function dataProviderForPositiveSemidefiniteMatrix(): array
    {
        return [
            [
                [
                    [0, 0],
                    [0, 0],
                ],
            ],
            [
                [
                    [1, 0],
                    [0, 1],
                ],
            ],
            [
                [
                    [1, 0],
                    [0, 2],
                ],
            ],
            [
                [
                    [1, 1],
                    [1, 1],
                ],
            ],
            [
                [
                    [2, -1],
                    [-1, 2],
                ],
            ],
            [
                [
                    [0, 0, 0],
                    [0, 3, 0],
                    [0, 0, 3],
                ],
            ],
            [
                [
                    [2, -1, -1],
                    [-1, 2, -1],
                    [-1, -1, 2],
                ],
            ],
            [
                [
                    [2, -1, 0],
                    [-1, 2, -1],
                    [0, -1, 2],
                ],
            ],
            [
                [
                    [2, -1, 1],
                    [-1, 2, -1],
                    [1, -1, 2],
                ],
            ],
            [
                [
                    [2, -1, 2],
                    [-1, 2, -1],
                    [2, -1, 2],
                ],
            ],
            [
                [
                    [9, -3, 3, 9],
                    [-3, 17, -1, -7],
                    [3, -1, 17, 15],
                    [9, -7, 15, 44],
                ],
            ],
        ];
    }

    /**
     * @testCase Axiom: Positive definiteness A is PD ⇒ A is PSD
     * If A is positive definite, then A is also positive semidefinite.
     * @dataProvider dataProviderForPositiveDefiniteMatrix
     * @param        array $A
     */
    public function testPositiveDefiniteIsAlsoPositiveSemidefinite(array $A)
    {
        $A  = MatrixFactory::create($A);

        $this->assertTrue($A->isPositiveDefinite());
        $this->assertTrue($A->isPositiveSemidefinite());
    }

    /**
     * @testCase Axiom: Negative definiteness A is ND ⇒ A is NSD
     * If A is negative definite, then A is also negative semidefinite.
     * @dataProvider dataProviderForNegativeDefiniteMatrix
     * @param        array $A
     */
    public function testNegativeDefiniteIsAlsoNegativeSemidefinite(array $A)
    {
        $A  = MatrixFactory::create($A);

        $this->assertTrue($A->isNegativeDefinite());
        $this->assertTrue($A->isNegativeSemidefinite());
    }

    public function dataProviderForNegativeDefiniteMatrix(): array
    {
        return [
            [
                [
                    [-1, 1],
                    [1, -2],
                ],
            ],
            [
                [
                    [-3, 0, 0],
                    [0, -2, 0],
                    [0, 0, -1],
                ],
            ],
        ];
    }

    /**
     * @testCase Axiom: Positive definiteness A is PD ⇔ A⁻¹ is PD
     * If A is positive definite, then A⁻¹ is positive definite.
     * @dataProvider dataProviderForPositiveDefiniteMatrix
     * @param        array $A
     */
    public function testPositiveDefiniteInverseIsPositiveDefinite(array $A)
    {
        $A   = MatrixFactory::create($A);
        $A⁻¹ = $A->inverse();

        // Floating point adjustment
        $A⁻¹ = $A⁻¹->map(function ($x) {
            return round($x, 7);
        });

        $this->assertTrue($A->isPositiveDefinite());
        $this->assertTrue($A⁻¹->isPositiveDefinite());
    }

    /**
     * @testCase Axiom: Negative definiteness A is ND ⇔ A⁻¹ is ND
     * If A is negative definite, then A⁻¹ is negative definite.
     * @dataProvider dataProviderForNegativeDefiniteMatrix
     * @param        array $A
     */
    public function testNegativeDefiniteInverseIsNegativeDefinite(array $A)
    {
        $A   = MatrixFactory::create($A);
        $A⁻¹ = $A->inverse();

        // Floating point adjustment
        $A⁻¹ = $A⁻¹->map(function ($x) {
            return round($x, 7);
        });

        $this->assertTrue($A->isNegativeDefinite());
        $this->assertTrue($A⁻¹->isNegativeDefinite());
    }

    /**
     * @testCase Axiom: Positive definiteness A is PD and r > 0 ⇒ rA is PD
     * If A is positive definite and r > 0, then rA is positive definite.
     * @dataProvider dataProviderForPositiveDefiniteMatrix
     * @param        array $A
     */
    public function testPositiveDefiniteThenScalarMultiplyWithPositiveNumberIsPositiveDefinite(array $A)
    {
        $A = MatrixFactory::create($A);
        $this->assertTrue($A->isPositiveDefinite());

        foreach (range(1, 10) as $r) {
            $rA = $A->scalarMultiply($r);
            $this->assertTrue($rA->isPositiveDefinite());
        }
        $A⁻¹ = $A->inverse();
    }

    /**
     * @testCase Axiom: Positive definiteness A and B are PD ⇒ A + B is PD
     * If A and B are positive definite then A + B is positive definite.
     * @dataProvider dataProviderForPositiveDefiniteMatrix
     * @param        array $A
     */
    public function testPositiveDefiniteAPlusAIsPositiveDefinite(array $M)
    {
        $A = MatrixFactory::create($M);
        $B = MatrixFactory::create($M);
        $A＋B = $A->add($B);

        $this->assertTrue($A->isPositiveDefinite());
        $this->assertTrue($B->isPositiveDefinite());
        $this->assertTrue($A＋B->isPositiveDefinite());
    }

    /**
     * @testCase Axiom: Positive definiteness A and B are PD ⇒ A + B is PD
     * If A and B are positive definite then A + B is positive definite.
     * @dataProvider dataProviderForTwoPositiveDefiniteMatrixes
     * @param  array $A
     * @param  array $B
     */
    public function testPositiveDefiniteAPlusBIsPositiveDefinite(array $A, array $B)
    {
        $A = MatrixFactory::create($A);
        $B = MatrixFactory::create($B);
        $A＋B = $A->add($B);

        $this->assertTrue($A->isPositiveDefinite());
        $this->assertTrue($B->isPositiveDefinite());
        $this->assertTrue($A＋B->isPositiveDefinite());
    }

    /**
     * @testCase Axiom: Positive definiteness A and B are PD ⇒ ABA is PD
     * If A and B are positive definite then ABA is positive definite.
     * @dataProvider dataProviderForPositiveDefiniteMatrix
     * @param  array $A
     */
    public function testPositiveDefiniteAAAIsPositiveDefinite(array $M)
    {
        $A = MatrixFactory::create($M);
        $B = MatrixFactory::create($M);
        $ABA = $A->multiply($B)->multiply($A);

        $this->assertTrue($A->isPositiveDefinite());
        $this->assertTrue($B->isPositiveDefinite());
        $this->assertTrue($ABA->isPositiveDefinite());
    }

    /**
     * @testCase Axiom: Positive definiteness A and B are PD ⇒ ABA is PD
     * If A and B are positive definite then ABA is positive definite.
     * @dataProvider dataProviderForTwoPositiveDefiniteMatrixes
     * @param  array $A
     * @param  array $B
     */
    public function testPositiveDefiniteABAIsPositiveDefinite(array $A, array $B)
    {
        $A = MatrixFactory::create($A);
        $B = MatrixFactory::create($B);
        $ABA = $A->multiply($B)->multiply($A);

        $this->assertTrue($A->isPositiveDefinite());
        $this->assertTrue($B->isPositiveDefinite());
        $this->assertTrue($ABA->isPositiveDefinite());
    }

    /**
     * @testCase Axiom: Positive definiteness A and B are PD ⇒ BAB is PD
     * If A and B are positive definite then BAB is positive definite.
     * @dataProvider dataProviderForTwoPositiveDefiniteMatrixes
     * @param  array $A
     * @param  array $B
     */
    public function testPositiveDefiniteBABIsPositiveDefinite(array $A, array $B)
    {
        $A = MatrixFactory::create($A);
        $B = MatrixFactory::create($B);
        $BAB = $B->multiply($A)->multiply($B);

        $this->assertTrue($A->isPositiveDefinite());
        $this->assertTrue($B->isPositiveDefinite());
        $this->assertTrue($BAB->isPositiveDefinite());
    }

    public function dataProviderForTwoPositiveDefiniteMatrixes(): array
    {
        return [
            [
                [
                    [2, -1],
                    [-1, 2],
                ],
                [
                    [1, -1],
                    [-1, 4],
                ],
            ],
            [
                [
                    [5, 2],
                    [2, 3],
                ],
                [
                    [6, 4],
                    [4, 5],
                ],
            ],
            [
                [
                    [12, -12],
                    [-12, 96],
                ],
                [
                    [6, 4],
                    [4, 5],
                ],
            ],
            [
                [
                    [2, -1, 0],
                    [-1, 2, -1],
                    [0, -1, 2],
                ],
                [
                    [2, -1, 1],
                    [-1, 2, -1],
                    [1, -1, 2],
                ],
            ],
            [
                [
                    [1, 0, 0],
                    [0, 3, 0],
                    [0, 0, 2],
                ],
                [
                    [3, -2, 0],
                    [-2, 2, 0],
                    [0, 0, 2],
                ],
            ],
            [
                [
                    [4, 1, -1],
                    [1, 2, 1],
                    [-1, 1, 2],
                ],
                [
                    [3, -2, 0],
                    [-2, 2, 0],
                    [0, 0, 2],
                ],
            ],
            [
                [
                    [14, 4, 9],
                    [4, 14, -7],
                    [9, -7, 14],
                ],
                [
                    [13, 0, -3],
                    [0, 9, 9],
                    [-3, 9, 10],
                ],
            ],
            [
                [
                    [14, -7, -13],
                    [-7, 6, 5],
                    [-13, 5, 14],
                ],
                [
                    [13, 0, -3],
                    [0, 9, 9],
                    [-3, 9, 10],
                ],
            ],
        ];
    }

    /**
     * @testCase Axiom: Zero matrix is lower triangular
     */
    public function testZeroMatrixIsLowerTriangular()
    {
        foreach (range(1, 20) as $m) {
            $L = MatrixFactory::zero($m, $m);
            $this->assertTrue($L->isLowerTriangular());
        }
    }

    /**
     * @testCase Axiom: Zero matrix is upper triangular
     */
    public function testZeroMatrixIsUpperTriangular()
    {
        foreach (range(1, 20) as $m) {
            $L = MatrixFactory::zero($m, $m);
            $this->assertTrue($L->isUpperTriangular());
        }
    }

    /**
     * @testCase Axiom: Zero matrix is diagonal
     */
    public function testZeroMatrixIsDiagonal()
    {
        foreach (range(1, 20) as $m) {
            $L = MatrixFactory::zero($m, $m);
            $this->assertTrue($L->isDiagonal());
        }
    }

    /**
     * @testCase Axiom: Lᵀ is upper triangular
     * Transpose of a lower triangular matrix is upper triagular
     * @dataProvider dataProviderForLowerTriangularMatrix
     * @param array $L
     */
    public function testTransposeOfLowerTriangularMatrixIsUpperTriangular(array $L)
    {
        $L  = MatrixFactory::create($L);
        $Lᵀ = $L->Transpose();

        $this->assertTrue($L->isLowerTriangular());
        $this->assertTrue($Lᵀ->isUpperTriangular());
    }

    /**
     * @testCase Axiom: Uᵀ is lower triangular
     * Transpose of an upper triangular matrix is lower triagular
     * @dataProvider dataProviderForUpperTriangularMatrix
     * @param array $U
     */
    public function testTransposeOfUpperTriangularMatrixIsLowerTriangular(array $U)
    {
        $U  = MatrixFactory::create($U);
        $Uᵀ = $U->Transpose();

        $this->assertTrue($U->isUpperTriangular());
        $this->assertTrue($Uᵀ->isLowerTriangular());
    }

    /**
     * @testCase Axiom: LL is lower triangular
     * Product of two lower triangular matrices is lower triangular
     * @dataProvider dataProviderForLowerTriangularMatrix
     * @param array $L
     */
    public function testProductOfTwoLowerTriangularMatricesIsLowerTriangular(array $L)
    {
        $L  = MatrixFactory::create($L);
        $LL = $L->multiply($L);

        $this->assertTrue($L->isLowerTriangular());
        $this->assertTrue($LL->isLowerTriangular());
    }

    /**
     * @testCase Axiom: UU is upper triangular
     * Product of two upper triangular matrices is upper triangular
     * @dataProvider dataProviderForUpperTriangularMatrix
     * @param array $U
     */
    public function testProductOfTwoUpperTriangularMatricesIsUpperTriangular(array $U)
    {
        $U  = MatrixFactory::create($U);
        $UU = $U->multiply($U);

        $this->assertTrue($U->isUpperTriangular());
        $this->assertTrue($UU->isUpperTriangular());
    }

    /**
     * @testCase Axiom: L + L is lower triangular
     * Sum of two lower triangular matrices is lower triangular
     * @dataProvider dataProviderForLowerTriangularMatrix
     * @param array $L
     */
    public function testSumOfTwoLowerTriangularMatricesIsLowerTriangular(array $L)
    {
        $L    = MatrixFactory::create($L);
        $L＋L = $L->add($L);

        $this->assertTrue($L->isLowerTriangular());
        $this->assertTrue($L＋L->isLowerTriangular());
    }

    /**
     * @testCase Axiom: U + U is upper triangular
     * Sum of two upper triangular matrices is upper triangular
     * @dataProvider dataProviderForUpperTriangularMatrix
     * @param array $U
     */
    public function testSumOfTwoUpperTriangularMatricesIsUpperTriangular(array $U)
    {
        $U    = MatrixFactory::create($U);
        $U＋U = $U->add($U);

        $this->assertTrue($U->isUpperTriangular());
        $this->assertTrue($U＋U->isUpperTriangular());
    }

    /**
     * @testCase Axiom: L⁻¹ is lower triangular (If L is invertible)
     * The inverse of an invertible lower triangular matrix is lower triangular
     * @dataProvider dataProviderForLowerTriangularMatrix
     * @param array $L
     */
    public function testInverseOfInvertibleLowerTriangularMatrixIsLowerTriangular(array $L)
    {
        $L = MatrixFactory::create($L);
        $this->assertTrue($L->isLowerTriangular());

        if ($L->isInvertible()) {
            $L⁻¹ = $L->inverse();
            $this->assertTrue($L⁻¹->isLowerTriangular());
        }
    }

    /**
     * @testCase Axiom: U⁻¹ is upper triangular (If U is invertible)
     * The inverse of an invertible upper triangular matrix is upper triangular
     * @dataProvider dataProviderForUpperTriangularMatrix
     * @param array $U
     */
    public function testInverseOfInvertibleUpperTriangularMatrixIsUpperTriangular(array $U)
    {
        $U = MatrixFactory::create($U);
        $this->assertTrue($U->isUpperTriangular());

        if ($U->isInvertible()) {
            $U⁻¹ = $U->inverse();
            $this->assertTrue($U⁻¹->isUpperTriangular());
        }
    }

    /**
     * @testCase Axiom: kL is lower triangular
     * Product of a lower triangular matrix by a constant is lower triangular
     * @dataProvider dataProviderForLowerTriangularMatrix
     * @param array $L
     */
    public function testProductOfLowerTriangularMatrixByConstantIsLowerTriangular(array $L)
    {
        $L = MatrixFactory::create($L);
        $this->assertTrue($L->isLowerTriangular());

        foreach (range(1, 10) as $k) {
            $kL = $L->scalarMultiply($k);
            $this->assertTrue($kL->isLowerTriangular());
        }
    }

    /**
     * @testCase Axiom: kU is upper triangular
     * Product of a upper triangular matrix by a constant is upper triangular
     * @dataProvider dataProviderForUpperTriangularMatrix
     * @param array $U
     */
    public function testProductOfUpperTriangularMatrixByConstantIsUpperTriangular(array $U)
    {
        $U = MatrixFactory::create($U);
        $this->assertTrue($U->isUpperTriangular());

        foreach (range(1, 10) as $k) {
            $kU = $U->scalarMultiply($k);
            $this->assertTrue($kU->isUpperTriangular());
        }
    }

    /**
     * @testCase Axiom: L is invertible iff diagonal is all non zero
     * Lower triangular matrix is invertible if and only if its diagonal entries are all non zero
     * @dataProvider dataProviderForLowerTriangularMatrix
     * @param array $L
     */
    public function testLowerTriangularMatrixIsInvertibleIfAndOnlyIfDigaonalEntriesAreAllNonZero(array $L)
    {
        $L = MatrixFactory::create($L);
        $this->assertTrue($L->isLowerTriangular());

        $zeros = array_filter(
            $L->getDiagonalElements(),
            function ($x) {
                return $x == 0;
            }
        );

        if (count($zeros) == 0) {
            $this->assertTrue($L->isInvertible());
        } else {
            $this->assertFalse($L->isInvertible());
        }
    }

    /**
     * @testCase Axiom: U is invertible iff diagonal is all non zero
     * Upper triangular matrix is invertible if and only if its diagonal entries are all non zero
     * @dataProvider dataProviderForUpperTriangularMatrix
     * @param array $U
     */
    public function testUpperTriangularMatrixIsInvertibleIfAndOnlyIfDigaonalEntriesAreAllNonZero(array $U)
    {
        $U = MatrixFactory::create($U);
        $this->assertTrue($U->isUpperTriangular());

        $zeros = array_filter(
            $U->getDiagonalElements(),
            function ($x) {
                return $x == 0;
            }
        );

        if (count($zeros) == 0) {
            $this->assertTrue($U->isInvertible());
        } else {
            $this->assertFalse($U->isInvertible());
        }
    }

    /**
     * @testCase Axiom: Dᵀ is diagonal
     * Transpose of a diagonal matrix is diagonal
     * @dataProvider dataProviderForDiagonalMatrix
     * @param array $D
     */
    public function testTransposeOfDiagonalMatrixIsDiagonal(array $D)
    {
        $D  = MatrixFactory::create($D);
        $Dᵀ = $D->Transpose();

        $this->assertTrue($D->isDiagonal());
        $this->assertTrue($Dᵀ->isDiagonal());
    }

    /**
     * @testCase Axiom: DD is diagonal
     * Product of two diagonal matrices is diagonal
     * @dataProvider dataProviderForDiagonalMatrix
     * @param array $D
     */
    public function testProductOfTwoDiagonalMatricesIsDiagonal(array $D)
    {
        $D  = MatrixFactory::create($D);
        $DD = $D->multiply($D);

        $this->assertTrue($D->isDiagonal());
        $this->assertTrue($DD->isDiagonal());
    }

    /**
     * @testCase Axiom: D + D is diagonal
     * Sum of two diagonal matrices is diagonal
     * @dataProvider dataProviderForDiagonalMatrix
     * @param array $D
     */
    public function testSumOfTwoDiagonalMatricesIsDiagonal(array $D)
    {
        $D    = MatrixFactory::create($D);
        $D＋D = $D->add($D);

        $this->assertTrue($D->isDiagonal());
        $this->assertTrue($D＋D->isDiagonal());
    }

    /**
     * @testCase Axiom: D⁻¹ is diagonal (If D is invertible)
     * The inverse of an invertible diagonal matrix is diagonal
     * @dataProvider dataProviderForDiagonalMatrix
     * @param array $D
     */
    public function testInverseOfInvertibleDiagonalMatrixIsDiagonal(array $D)
    {
        $D = MatrixFactory::create($D);
        $this->assertTrue($D->isDiagonal());

        if ($D->isInvertible()) {
            $D⁻¹ = $D->inverse();
            $this->assertTrue($D⁻¹->isDiagonal());
        }
    }

    /**
     * @testCase Axiom: kD is Diagonal
     * Product of a diagonal matrix by a constant is diagonal
     * @dataProvider dataProviderForDiagonalMatrix
     * @param array $D
     */
    public function testProductOfDiagonalMatrixByConstantIsDiagonal(array $D)
    {
        $D = MatrixFactory::create($D);
        $this->assertTrue($D->isDiagonal());

        foreach (range(1, 10) as $k) {
            $kD = $D->scalarMultiply($k);
            $this->assertTrue($kD->isDiagonal());
        }
    }

    /**
     * @testCase Axiom: D is invertible iff diagonal is all non zero
     * Diagonal matrix is invertible if and only if its diagonal entries are all non zero
     * @dataProvider dataProviderForDiagonalMatrix
     * @param array $D
     */
    public function testDiagonalMatrixIsInvertibleIfAndOnlyIfDigaonalEntriesAreAllNonZero(array $D)
    {
        $D = MatrixFactory::create($D);
        $this->assertTrue($D->isDiagonal());

        $zeros = array_filter(
            $D->getDiagonalElements(),
            function ($x) {
                return $x == 0;
            }
        );

        if (count($zeros) == 0) {
            $this->assertTrue($D->isInvertible());
        } else {
            $this->assertFalse($D->isInvertible());
        }
    }

    public function dataProviderForLowerTriangularMatrix(): array
    {
        return [
            [
                [
                    [1],
                ],
            ],
            [
                [
                    [0],
                ],
            ],
            [
                [
                    [1, 0],
                    [1, 1],
                ],
            ],
            [
                [
                    [1, 0, 0],
                    [1, 1, 0],
                    [1, 1, 1],
                ],
            ],
            [
                [
                    [1, 0, 0],
                    [2, 3, 0],
                    [4, 5, 6],
                ],
            ],
            [
                [
                    [1, 0, 0],
                    [2, 0, 0],
                    [4, 5, 6],
                ],
            ],
            [
                [
                    [1, 0, 0, 0],
                    [1, 1, 0, 0],
                    [1, 1, 1, 0],
                    [1, 1, 1, 1],
                ],
            ],
            [
                [
                    [5, 0, 0, 0],
                    [-6, 1, 0, 0],
                    [4, 6, 8, 0],
                    [6, 7, 7, -1],
                ],
            ],
            [
                [
                    [5, 0, 0, 0],
                    [-6, 1, 0, 0],
                    [4, 6, 0, 0],
                    [6, 7, 7, -1],
                ],
            ],
            [
                [
                    [1, 0, 0, 0, 0, 0],
                    [1, 2, 0, 0, 0, 0],
                    [1, 0, 3, 0, 0, 0],
                    [1, 2, 0, 4, 0, 0],
                    [1, 0, 0, 0, 5, 0],
                    [1, 2, 3, 0, 0, 6],
                ],
            ],
        ];
    }

    public function dataProviderForUpperTriangularMatrix(): array
    {
        return [
            [
                [
                    [1],
                ],
            ],
            [
                [
                    [0],
                ],
            ],
            [
                [
                    [1, 1],
                    [0, 1],
                ]
            ],
            [
                [
                    [1, 2],
                    [0, 4],
                ],
            ],
            [
                [
                    [1, 2, 3],
                    [0, 4, 5],
                    [0, 0, 6],
                ],
            ],
            [
                [
                    [6, 5, 4],
                    [0, 8, 8],
                    [0, 0, 9],
                ],
            ],
            [
                [
                    [6, 5, 4],
                    [0, 8, 8],
                    [0, 0, 0],
                ],
            ],
            [
                [
                    [1, 2, 3, 4],
                    [0, 4, 5, 6],
                    [0, 0, 6, 7],
                    [0, 0, 0, 8],
                ],
            ],
            [
                [
                    [-1, 0, 0, 0, 0, 0],
                    [0, 2, 0, 0, 0, 0],
                    [0, 0, 3, 0, 0, 0],
                    [0, 0, 0, 4, 0, 0],
                    [0, 0, 0, 0, 5, 0],
                    [0, 0, 0, 0, 0, 6],
                ],
            ],
            [
                [
                    [1, 0, 0, 0, 0, 0],
                    [0, 2, 0, 0, 0, 0],
                    [0, 0, 3, 0, 0, 0],
                    [0, 0, 0, 4, 0, 0],
                    [0, 0, 0, 0, 5, 0],
                    [0, 0, 0, 0, 0, 6],
                ],
            ],
            [
                [
                    [1, 0, 0, 0, 0, 0],
                    [0, 2, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 5, 0],
                    [0, 0, 0, 0, 0, 6],
                ],
            ],
        ];
    }

    public function dataProviderForDiagonalMatrix(): array
    {
        return [
            [
                [
                    [0],
                ],
            ],
            [
                [
                    [1],
                ],
            ],
            [
                [
                    [1, 0],
                    [0, 1],
                ],
            ],
            [
                [
                    [-5, 0],
                    [0, 3],
                ],
            ],
            [
                [
                    [1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1],
                ],
            ],
            [
                [
                    [1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1],
                ],
            ],
            [
                [
                    [1, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0],
                    [0, 0, 1, 0, 0],
                    [0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 1],
                ],
            ],
            [
                [
                    [1, 0, 0, 0, 0],
                    [0, 2, 0, 0, 0],
                    [0, 0, 3, 0, 0],
                    [0, 0, 0, 4, 0],
                    [0, 0, 0, 0, -6],
                ],
            ],
            [
                [
                    [-1, 0, 0, 0, 0, 0],
                    [0, 2, 0, 0, 0, 0],
                    [0, 0, 3, 0, 0, 0],
                    [0, 0, 0, 4, 0, 0],
                    [0, 0, 0, 0, 5, 0],
                    [0, 0, 0, 0, 0, 6],
                ],
            ],
            [
                [
                    [1, 0, 0, 0, 0, 0],
                    [0, 2, 0, 0, 0, 0],
                    [0, 0, 3, 0, 0, 0],
                    [0, 0, 0, 4, 0, 0],
                    [0, 0, 0, 0, 5, 0],
                    [0, 0, 0, 0, 0, 6],
                ],
            ],
            [
                [
                    [1, 0, 0, 0, 0, 0],
                    [0, 2, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 4, 0, 0],
                    [0, 0, 0, 0, 5, 0],
                    [0, 0, 0, 0, 0, 0],
                ],
            ],
        ];
    }

    /**
     * @testCase Axiom: Reduced row echelon form is upper triangular
     * @dataProvider dataProviderForOneSquareMatrix
     * @param array $A
     */
    public function testReducedRowEchelonFormIsUpperTriangular(array $A)
    {
        $A    = MatrixFactory::create($A);
        $rref = $A->rref();

        $this->assertTrue($rref->isUpperTriangular());
    }
}
