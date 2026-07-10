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
        # 원소 기호 읽기 (대문자로 시작, 소문자 이어질 수 있음)
        element = formula[i]
        i += 1
        if i < len(formula) and formula[i].islower():
            element += formula[i]
            i += 1
        
        # 숫자 읽기 (개수)
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
    
    print(f"\n화학식: {formula}")
    print("-" * 30)
    
    for element, count in elements.items():
        if element not in atomic_weights:
            print(f"⚠️ '{element}' 원소의 원자량 정보가 없습니다.")
            return None
        
        weight = atomic_weights[element] * count
        total_weight += weight
        print(f"{element}: {count}개 x {atomic_weights[element]} = {weight:.3f}")
    
    print("-" * 30)
    print(f"총 분자량: {total_weight:.3f} g/mol")
    return total_weight


# 프로그램 시작
print("=== 분자량 계산기 ===")
print("예시: H2O (물), C6H12O6 (포도당), NaCl (소금)")
print()

while True:
    user_formula = input("화학식을 입력하세요 (종료하려면 'q' 입력): ")
    
    if user_formula.lower() == 'q':
        print("프로그램을 종료합니다.")
        break
    
    calculate_molecular_weight(user_formula)
    print()