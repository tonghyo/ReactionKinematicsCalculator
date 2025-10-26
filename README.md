# Fusion Reaction Simulation

이 프로젝트는 핵융합 반응 시뮬레이션을 위한 C++ 프로그램입니다.

## 파일 구조

- `FusionReaction.h` - 클래스 선언 및 헤더 파일
- `FusionReaction_Setup.cpp` - 설정 함수들 (SetBeamParameters, SetTargetParameters, AddProduct)
- `FusionReaction_MassHist.cpp` - 질량 파일 읽기 및 히스토그램 초기화
- `FusionReaction_Kinematics.cpp` - 운동학 계산 함수들
- `FusionReaction_Analysis.cpp` - 분석 및 시뮬레이션 함수들
- `fusion_reaction.C` - 메인 실행 파일
- `Makefile` - 컴파일 설정
- `mass.dat` - 핵종 질량 데이터

## 컴파일 및 실행

### 요구사항
- ROOT (CERN의 데이터 분석 프레임워크)
- C++11 이상 지원 컴파일러

### 컴파일
```bash
make
```

### 실행
```bash
make run
```

### 정리
```bash
make clean
```

## 사용법

### 기본 설정
```cpp
FusionReaction reaction;

// 빔 설정 (에너지, A, Z)
reaction.SetBeamParameters(85.0, 17, 9); // 85 MeV 17F

// 타겟 설정 (A, Z)
reaction.SetTargetParameters(28, 14); // 28Si

// 생성물 추가 (A, Z, 이름)
reaction.AddProduct(42, 22, "42V");   // 42V
reaction.AddProduct(1, 0, "n1");      // neutron 1
reaction.AddProduct(1, 0, "n2");      // neutron 2
reaction.AddProduct(1, 0, "n3");      // neutron 3
```

### 시뮬레이션 실행
```cpp
// 질량 파일 읽기
reaction.ReadMassFile("mass.dat");

// 히스토그램 초기화
reaction.InitializeHistograms();

// 시뮬레이션 실행 (이벤트 수, 상세 출력 여부)
reaction.RunSimulation(10000, true);

// 결과 저장
reaction.SaveResults("fusion_results.root");
```

## 출력 파일

- `fusion_results.root` - ROOT 형식의 히스토그램 파일
  - 빔 에너지 분포
  - 생성물 각도 및 에너지 분포
  - Lab frame vs CM frame 비교
  - 에너지 보존 검증

## 주요 기능

1. **다체 반응 시뮬레이션**: 임의의 수의 생성물을 가진 핵반응
2. **상대론적 운동학**: 정확한 Lorentz 변환
3. **실험적 해상도**: 각도 및 에너지 측정 불확실성 반영
4. **에너지 보존 검증**: 시뮬레이션 정확성 확인
5. **다양한 히스토그램**: CM frame과 Lab frame 비교 분석
