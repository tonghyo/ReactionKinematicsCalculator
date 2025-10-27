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

## Usage notes (parameter file)

이 프로젝트는 이제 간단한 key=value 형식의 파라미터 파일을 통해 시뮬레이션 파라미터를 로드할 수 있습니다.
기본 파일명은 `params.txt`이며, 프로젝트 루트에 두거나 컴파일된 바이너리를 실행할 때 첫 번째 인자로 파일 경로를 넘기면 됩니다.

예제 파일 `params_example.txt`가 리포지토리에 포함되어 있습니다. 주요 키 예시는 아래와 같습니다:

- `beam` = Energy,A,Z
- `target` = A,Z
- `experimental` = E_loss,E_strag,E_beam_re,tar_res,th_res_deg
- `products` = A,Z,label;A,Z,label;...
- `excited_energies`, `excited_branching` = comma-separated lists (같은 길이여야 함)
- `mass_file`, `n_events`, `output_file`, `verbose_events`, `no_draw` 등

ROOT에서 매크로로 호출하면 기본적으로 `params.txt`를 참조합니다. 컴파일된 실행파일을 사용할 때는 파라미터 파일 경로를 넘기세요:

  ./fusion_reaction params_example.txt


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
