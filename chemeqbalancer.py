import re

import numpy as np
import streamlit as st
from sympy import Matrix, lcm

ELEMENTS = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "K": 19,
    "Ca": 20,
}


def parse_elements(equation):
    elements = set(re.findall(r"[A-Z][a-z]?", equation))
    for el in elements:
        if el not in ELEMENTS:
            return False, el
    return True, None


def parse_formula(formula):
    tokens = re.findall(r"([A-Z][a-z]?)(\d*)", formula)
    counts = {}
    for el, cnt in tokens:
        if el not in ELEMENTS:
            raise ValueError(f"지원하지 않는 원소입니다: {el}")
        counts[el] = counts.get(el, 0) + int(cnt) if cnt else counts.get(el, 0) + 1
    return counts


def parse_side(side):
    return [s.strip() for s in side.split("+")]


def get_element_list(compounds):
    elements = set()
    for cmpd in compounds:
        elements.update(parse_formula(cmpd).keys())
    return sorted(elements, key=lambda x: ELEMENTS[x])


def build_matrix(lhs, rhs, elements):
    matrix = []
    for el in elements:
        row = []
        for cmpd in lhs:
            row.append(parse_formula(cmpd).get(el, 0))
        for cmpd in rhs:
            row.append(-parse_formula(cmpd).get(el, 0))
        matrix.append(row)
    return np.array(matrix)


def balance(equation):
    lhs, rhs = equation.split("->")
    lhs = parse_side(lhs)
    rhs = parse_side(rhs)
    elements = get_element_list(lhs + rhs)
    mat = build_matrix(lhs, rhs, elements)
    M = Matrix(mat)
    nullsp = M.nullspace()
    if not nullsp:
        return None
    vec = nullsp[0]
    lcm_val = lcm([r.q for r in vec])
    coeffs = [abs(int(r * lcm_val)) for r in vec]
    left = " + ".join(
        f"{coeff if coeff > 1 else ''}{cmpd}"
        for coeff, cmpd in zip(coeffs[: len(lhs)], lhs)
    )
    right = " + ".join(
        f"{coeff if coeff > 1 else ''}{cmpd}"
        for coeff, cmpd in zip(coeffs[len(lhs) :], rhs)
    )
    return f"{left} -> {right}"


def main():
    st.title("화학 반응식 균형 맞추기")
    st.write("H~Ca까지의 원소만 지원합니다.")
    eq = st.text_input("화학 반응식을 입력하세요 (예: NH3 + O2 -> NO + H2O)")
    if eq:
        valid, el = parse_elements(eq)
        if not valid:
            st.error(f"지원하지 않는 원소입니다: {el}")
            return
        try:
            result = balance(eq.replace(" ", ""))
            if result:
                st.success(f"균형 맞춘 결과: {result}")
            else:
                st.error("균형을 맞출 수 없습니다.")
        except Exception as e:
            st.error(f"오류: {e}")


if __name__ == "__main__":
    main()
