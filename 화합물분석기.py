from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

def analyze_compound(smiles):
    """SMILES 문자열로 화합물을 분석하는 함수"""
    
    # SMILES를 분자 객체로 변환
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        print("⚠️ 올바른 SMILES 표기법이 아닙니다.")
        return
    
    print(f"\n입력한 SMILES: {smiles}")
    print("-" * 40)
    
    # 기본 정보 계산
    mol_weight = Descriptors.MolWt(mol)  # 분자량
    num_atoms = mol.GetNumAtoms()  # 원자 개수 (수소 제외)
    num_h_donors = Descriptors.NumHDonors(mol)  # 수소 결합 공여체
    num_h_acceptors = Descriptors.NumHAcceptors(mol)  # 수소 결합 수용체
    logp = Descriptors.MolLogP(mol)  # 지용성/수용성 지표
    
    print(f"분자량: {mol_weight:.2f} g/mol")
    print(f"원자 개수 (수소 제외): {num_atoms}개")
    print(f"수소 결합 공여체 수: {num_h_donors}")
    print(f"수소 결합 수용체 수: {num_h_acceptors}")
    print(f"LogP (지용성 지표): {logp:.2f}")
    print("-" * 40)
    
    # 신약 후보로서의 간단한 평가 (Lipinski's Rule of 5)
    print("💊 신약 후보 가능성 평가 (Lipinski's Rule of 5):")
    rules_passed = 0
    
    if mol_weight <= 500:
        print("  ✅ 분자량 500 이하")
        rules_passed += 1
    else:
        print("  ❌ 분자량 500 초과")
    
    if num_h_donors <= 5:
        print("  ✅ 수소 결합 공여체 5개 이하")
        rules_passed += 1
    else:
        print("  ❌ 수소 결합 공여체 5개 초과")
    
    if num_h_acceptors <= 10:
        print("  ✅ 수소 결합 수용체 10개 이하")
        rules_passed += 1
    else:
        print("  ❌ 수소 결합 수용체 10개 초과")
    
    if logp <= 5:
        print("  ✅ LogP 5 이하")
        rules_passed += 1
    else:
        print("  ❌ LogP 5 초과")
    
    print(f"\n{rules_passed}/4 조건 통과")
    if rules_passed >= 3:
        print("👍 경구 투여 가능한 약물이 될 가능성이 높습니다!")
    else:
        print("⚠️ 경구 투여용으로는 어려울 수 있습니다.")
    
    # 분자 구조 이미지로 저장
    img = Draw.MolToImage(mol, size=(400, 400))
    img.save("분자구조.png")
    print("\n📷 분자 구조 이미지가 '분자구조.png'로 저장되었습니다.")


# 프로그램 시작
print("=== 화합물 분석기 (SMILES 기반) ===")
print("예시:")
print("  아스피린: CC(=O)OC1=CC=CC=C1C(=O)O")
print("  카페인:   CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
print("  에탄올:   CCO")
print()

while True:
    smiles_input = input("SMILES를 입력하세요 (종료하려면 'q' 입력): ")
    
    if smiles_input.lower() == 'q':
        print("프로그램을 종료합니다.")
        break
    
    analyze_compound(smiles_input)
    print()