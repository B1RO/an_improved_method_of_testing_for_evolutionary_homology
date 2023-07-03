import {useState} from 'react'
import reactLogo from './assets/react.svg'
import './App.css'
import MatrixComponent from "./MatrixComponent";

enum AminoAcid {
    D = 'D',
    C = 'C',
    T = 'T',
    F = 'F',
    E = 'E',
    H = 'H',
    K = 'K',
    A = 'A',
    M = 'M',
    N = 'N',
    Y = 'Y',
    P = 'P',
    Q = 'Q',
    R = 'R',
    S = 'S',
    W = 'W',
    L = 'L',
    V = 'V',
    I = 'I',
    G = 'G',
}

enum NucleicAcid {
    Adenine = 'A',
    Thymine = 'T',
    Cytosine = 'C',
    Guanine = 'G',
}

type Codon = `${NucleicAcid}${NucleicAcid}${NucleicAcid}`;

type CodonTable = Record<Codon, AminoAcid | null>;

const codonTable: CodonTable = {
    AAA: AminoAcid.K,
    AAC: AminoAcid.N,
    AAG: AminoAcid.K,
    AAT: AminoAcid.N,
    ACA: AminoAcid.T,
    ACC: AminoAcid.T,
    ACG: AminoAcid.T,
    ACT: AminoAcid.T,
    AGA: AminoAcid.R,
    AGC: AminoAcid.S,
    AGG: AminoAcid.R,
    AGT: AminoAcid.S,
    ATA: AminoAcid.I,
    ATC: AminoAcid.I,
    ATG: AminoAcid.M,
    ATT: AminoAcid.I,
    CAA: AminoAcid.Q,
    CAC: AminoAcid.H,
    CAG: AminoAcid.Q,
    CAT: AminoAcid.H,
    CCA: AminoAcid.P,
    CCC: AminoAcid.P,
    CCG: AminoAcid.P,
    CCT: AminoAcid.P,
    CGA: AminoAcid.R,
    CGC: AminoAcid.R,
    CGG: AminoAcid.R,
    CGT: AminoAcid.R,
    CTA: AminoAcid.L,
    CTC: AminoAcid.L,
    CTG: AminoAcid.L,
    CTT: AminoAcid.L,
    GAA: AminoAcid.E,
    GAC: AminoAcid.D,
    GAG: AminoAcid.E,
    GAT: AminoAcid.D,
    GCA: AminoAcid.A,
    GCC: AminoAcid.A,
    GCG: AminoAcid.A,
    GCT: AminoAcid.A,
    GGA: AminoAcid.G,
    GGC: AminoAcid.G,
    GGG: AminoAcid.G,
    GGT: AminoAcid.G,
    GTA: AminoAcid.V,
    GTC: AminoAcid.V,
    GTG: AminoAcid.V,
    GTT: AminoAcid.V,
    TAA: null,
    TAC: AminoAcid.Y,
    TAG: null,
    TAT: AminoAcid.Y,
    TCA: AminoAcid.S,
    TCC: AminoAcid.S,
    TCG: AminoAcid.S,
    TCT: AminoAcid.S,
    TGA: null,
    TGC: AminoAcid.C,
    TGG: AminoAcid.W,
    TGT: AminoAcid.C,
    TTA: AminoAcid.L,
    TTC: AminoAcid.F,
    TTG: AminoAcid.L,
    TTT: AminoAcid.F,
};


function cartesianProduct<T>(arr1: T[], arr2: T[]): Array<[T, T]> {
    return arr1.flatMap((item1) => arr2.map((item2) => {
        const toReturn: [T, T] = [item1, item2]
        return toReturn;
    }));
}


function countDifferences([triplet1, triplet2]: Array<string>): number {
    return triplet1.split('').reduce((count, char, index) => count + (char !== triplet2[index] ? 1 : 0), 0);
}

function computeReplacementCount(aminoAcidA: string, aminoAcidB: string): number {
    let possibleSourceCodons = Object.entries(codonTable)
        .filter(([codon, aminoAcid]) => {
            return aminoAcid == aminoAcidA;
        })
        .map(([codon, aminoAcid]) => codon);
    let possibleTargetCodons = Object.entries(codonTable)
        .filter(([codon, aminoAcid]) => aminoAcid == aminoAcidB)
        .map(([codon, aminoAcid]) => codon);

    return Math.min(...cartesianProduct(possibleSourceCodons, possibleTargetCodons).map(countDifferences));
}

function generateReplacementMatrix() {
    let replacementMatrix: Array<Array<number>> = []
    for (let i = 0; i < 20; i++) {
        replacementMatrix[i] = [];
        for (let j = 0; j < 20; j++) {
            replacementMatrix[i][j] = computeReplacementCount(Object.keys(AminoAcid)[i], Object.keys(AminoAcid)[j]);
        }
    }
    return replacementMatrix
}

const getAllPossibleSubsequencesOfLength = (sequence: string, subsequenceLength: number) =>
    [...Array(sequence.length - subsequenceLength + 1).keys()].map(subsequenceStart => sequence.slice(subsequenceStart, subsequenceStart + subsequenceLength))

const computeReplacementCost = (costMatrix: Array<Array<number>>) => ([a, b]: [string, string]) => {
    return [...Array(a.length).keys()].map(index => {
        const aminoAcidIndexA = Object.keys(AminoAcid).findIndex(x => x == a[index])
        const aminoAcidIndexB = Object.keys(AminoAcid).findIndex(x => x == b[index])
        return costMatrix[aminoAcidIndexA][aminoAcidIndexB]
    }).reduce((a,b)=>a+b,0)
}

function factorial(n: number): number {
    if (n <= 1)
        return 1;
    else
        return n * factorial(n - 1)

}

function computeProbabilityOfOccurrence(replacementCosts: Array<number>, {p0, p1, p2, p3}: {
    p0: number,
    p1: number,
    p2: number,
    p3: number
}) {
    let LSE = replacementCosts.length;
    let a = replacementCosts.filter(x => x == 0).length
    let b = replacementCosts.filter(x => x == 1).length
    let c = replacementCosts.filter(x => x == 2).length
    let d = replacementCosts.filter(x => x == 3).length
    return Math.pow(p0 ,a)* Math.pow(p1,b) * Math.pow(p2,c) * Math.pow(p3,d) * factorial(LSE) / (factorial(a) * factorial(b) * factorial(c) * factorial(d))
}

function VMatrix(sequenceA : string, sequenceB : string)
{
    let frequencyMatrix: Array<Array<number>> = []
    for (let i = 0; i < 20; i++) {
        frequencyMatrix[i] = [];
        for (let j = 0; j < 20; j++) {
            let ithAminoAcid = Object.keys(AminoAcid)[i]
            let jthAminoAcid = Object.keys(AminoAcid)[j]
            let f_i = [...sequenceA].filter(char=>char==ithAminoAcid).length/sequenceA.length
            let f_j = [...sequenceB].filter(char=>char==jthAminoAcid).length/sequenceB.length
            frequencyMatrix[i][j]= f_i*f_j;
        }
    }
    return frequencyMatrix
}

function App() {
    const alphaChain = "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
    const betaChain = "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
    const LSE = 30
    const costMatrix = generateReplacementMatrix();
    let subsequencePairs: Array<[string, string]> = cartesianProduct<string>(getAllPossibleSubsequencesOfLength(alphaChain, 30), getAllPossibleSubsequencesOfLength(betaChain, 30))
    let subsequencePairsWithCost = subsequencePairs.map((pair) => ({
        pair: pair,
        cost: computeReplacementCost(costMatrix)(pair)
    }))

    console.log(subsequencePairsWithCost.filter(x=>x.cost<=34))
    return (
        <div>
            <MatrixComponent matrix={costMatrix}/>
        </div>
    )
}

export default App
