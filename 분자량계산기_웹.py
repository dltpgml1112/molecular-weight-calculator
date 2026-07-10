import streamlit as st

# 원소별 원자량 (단위: g/mol)
atomic_weights = {
    "H": 1.008,
    "C": 12.011,
    "O": 15.999,
    "N": 14.007,
    "Na": 22.990,
    "Cl": 35.453,
    "S": 32.06,
}

def parse_formula(formula):
    """화학식을 원소와 개수로 분리하는 함수"""
    result = {}
    i = 0
    while i < len(formula):
        element = formula[i]
        i += 1
        if i < len(formula) and formula[i].islower():
            element += formula[i]
            i += 1
        
        num = ""
        while i < len(formula) and formula[i].isdigit():
            num += formula[i]
            i += 1
        
        count = int(num) if num else 1
        
        if element in result:
            result[element] += count
        else:
            result[element] = count
    
    return result

def calculate_molecular_weight(formula):
    """분자량을 계산하는 함수"""
    elements = parse_formula(formula)
    total_weight = 0
    details = []
    
    for element, count in elements.items():
        if element not in atomic_weights:
            return None, None, element  # 모르는 원소
        
        weight = atomic_weights[element] * count
        total_weight += weight
        details.append((element, count, atomic_weights[element], weight))
    
    return total_weight, details, None


# ===== 웹 화면 구성 시작 =====

st.title("🧪 분자량 계산기")
st.write("화학식을 입력하면 분자량을 계산해드려요.")

# 예시 보여주기
st.info("예시: H2O (물), C6H12O6 (포도당), NaCl (소금)")

# 입력창
formula = st.text_input("화학식을 입력하세요:")

# 계산 버튼
if st.button("계산하기"):
    if formula:
        total, details, unknown_element = calculate_molecular_weight(formula)
        
        if unknown_element:
            st.error(f"⚠️ '{unknown_element}' 원소의 정보가 없습니다.")
        else:
            st.success(f"✅ 총 분자량: **{total:.3f} g/mol**")
            
            st.write("### 상세 계산")
            for element, count, weight_per, total_w in details:
                st.write(f"- {element}: {count}개 × {weight_per} = {total_w:.3f}")
    else:
        st.warning("화학식을 입력해주세요!")