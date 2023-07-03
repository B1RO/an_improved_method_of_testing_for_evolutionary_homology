import React from 'react';
import styled from 'styled-components';

interface MatrixProps {
    matrix: number[][];
}

const headers = ['A', 'C', 'E', 'F', 'G', 'H', 'I', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'Y'];

const Table = styled.table`
  border-collapse: collapse;
  width: 100%;
`;

const TableHeader = styled.th`
  padding: 10px;
  background-color: #f2f2f2;
`;

const TableData = styled.td`
  padding: 10px;
  text-align: center;
`;

export const MatrixComponent: React.FC<MatrixProps> = ({ matrix }) => {
    return (
        <Table>
            <thead>
            <tr>
                <TableHeader></TableHeader>
                {headers.map((header) => (
                    <TableHeader key={header}>{header}</TableHeader>
                ))}
            </tr>
            </thead>
            <tbody>
            {matrix.map((row, rowIndex) => (
                <tr key={rowIndex}>
                    <TableHeader>{headers[rowIndex]}</TableHeader>
                    {row.map((cell, cellIndex) => (
                        <TableData key={cellIndex}>{cell}</TableData>
                    ))}
                </tr>
            ))}
            </tbody>
        </Table>
    );
};

export default MatrixComponent;
