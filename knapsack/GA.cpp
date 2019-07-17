#include<iostream>
#include<algorithm>
#include<stdlib.h>
#include<time.h>
#include <windows.h>
#include <stdio.h>
#include<utility>
using namespace std;
#define CHECK_TIME_START __int64 freq, start, end; if (QueryPerformanceFrequency((_LARGE_INTEGER*)&freq)) {QueryPerformanceCounter((_LARGE_INTEGER*)&start);
#define CHECK_TIME_END(a,b) QueryPerformanceCounter((_LARGE_INTEGER*)&end); a=(float)((double)(end - start)/freq*1000); b=TRUE; } else b=FALSE;
#define CL 100 // chronomosome length 유전자 길이
#define PS 1000 // 한 세대의 인구
#define GN 5000 // 세대의 수(반복 수)
#define CAPACITY 3001// 최대 수용무게
#define CP1 0.3 // 교차point1
#define CP2 0.7 // 교차 point2
#define PC 0.9 // 교차가 발생하는 확률 90%
#define MP 0.05 // 변이율 1%
#define SMP 0.07 // 선택적 변이율 (SMP*100)% 를 랜덤하게 0으로 만든다.
#define EP 0.04 // 엘리티즘 비율, 상위 (EP*100)%를 다음세대에 포함시킨다.
#define FC 3 // 적응도 상수 (k > 1) 높을수록 선택압 증가

// 아이템을 담을 구조체
struct item
{
	short weight;
	short price;
};

// 각 세대의, 각 유전자에 대한 정보를 담을 구조체
struct P
{
	int chromo_num = 0;
	int best_profit = 0;
	int weight = 0;
	bool selective_mutate = false;
};

// chromosomes[세대][인구][유전자 길이]
// 세대 : 0 -> 현재 세대
//			  1 -> 다음 세대
int chromosomes[2][PS][CL];


int selective_mutated_count = 0; // 선택변이 발생카운터
int mutated_count = 0; // 변이 발생카운터
int fitness[PS]; // 적응도 배열
bool visited[PS+1] = { false }; // 다음 유전자 방문배열, 랜덤하게 삽입할때 중복 방지용
typedef pair<int, int> C; // 1. index, 2. fitness -> 적응도 정렬시 사용, in-place 유지를 위해 pair 사용
C cpyFitness[PS]; // pair C형 배열을 생성함으로써 적응도를 in-place하게 정렬 가능

void make_item(); // item list를 만든다
void read_item(item *items); // item list를 불러온다
void make_chromosomes(); // 초기 유전자를 생성한다 랜덤으로 0,1을 부여한다.
void print_chromosomes(); // 유전자를 출력한다.
P determine_fitness(item *items, int generation); // 적응도를 결정한다.
void print_fitness(int geneartion); // 적응도를 출력한다.
int rulet(int chromo); // 룰렛 함수
void cross(int generation, int chromo1, int chromo2); // chromo1과 chromo2의 교차 함수
void print_cross(int generation, int chromo1, int chromo2); // 교차의 결과를 출력한다.
void copy_chromosomes(int generation, int start, int chromo1, int chromo2); // 유전자를 다음세대로 복제한다.
void print_generation(int gen_start, int generation); // 세대 전체를 출력한다.
void mutation(int generation, int chromo1); // 변이함수.
void seletive_mutation(int generation, int chromo1); // 선택적 변이 함수
int elitism_select(int start); // 엘리티즘 함수
void copy_fitness(); // 적응도 배열 복사
void insertion_sort(C *data, int n); // 삽입정렬로 복사된 적응도 배열을 정렬
void init_visited(); // 방문배열을 모두 false로 초기화



int main()
{
	// item 배열
	int start = 0; // 유전자 반복 Count
	int gen_start = 0; // 세대 반복 Count
	P nextGen_Fitness; // 적응도 평가에 사용
	P bestChrom; // 최적의 유전자
	int bestGen = 0; // 최적의 유전자가 속한 세대
	int chromo1; // 선택 유전자 1
	int chromo2; // 선택 유전자 2 
	int generation = 0;
	item *items;
	items = (item *)malloc(sizeof(item)*CL);

	// 아이템 초기화
	read_item(items);

	// 유전자 생성
	make_chromosomes();

	// 초기 적응도 검사
	P init_gen = determine_fitness(items, generation);

	// 생성 사이클 시작
	while (gen_start < GN)
	{
		// 다음 세대 생성 시작
		for (start = 0; start < PS - 1; start += 2)
		{
			// EP 비율만큼의 유전자를
			// 엘리티즘으로 먼저 포함될 유전자로 선별
			if (start < EP*PS)
			{
				copy_fitness(); // 적응도를 복사 한뒤
				insertion_sort(cpyFitness, PS); // 적응도 정렬
				chromo1 = elitism_select(start + 100); // 적응도에 따라 선택
				chromo2 = elitism_select(start + 50000); // 인자값은 random seed
			}
			else // 엘리티즘 유전자 선별이 끝나면 룰렛으로 선별
			{
				chromo1 = rulet(start + gen_start+100); // 적응도 배열을 기반으로 한 룰렛함수
				chromo2 = rulet(start); // 인자값은 random seed
			}

			// 교차
			srand(time(NULL) + gen_start + start); // random seed
			if (rand() % 100 <= PC * 100)
			{
				cross(generation, chromo1, chromo2); // 선택된 chromo1 과 chromo2를 교차
			}
			
			// 다음 세대 시작 전, 다음 세대에 대한 변이
			// 가방의 무게를 넘는 유전자에 대해서는 선택적 변이를 수행
			if (nextGen_Fitness.selective_mutate)
			{
				seletive_mutation(generation, start); // 선택적 변이는 반전이 아니라 0으로 만든다.
				seletive_mutation(generation, start + 1); // start가 2번씩 뒤므로 2번씩
			}
			else
			{
				// 가방의 무게를 넘지 않는 경우 변이 수행
				mutation(generation, start);
				mutation(generation, start + 1);
			}

			// 다음세대에 복제
			copy_chromosomes(generation, start, chromo1, chromo2);
		}
		init_visited(); // 방문 배열 false로 초기화
		generation = !generation; // 다음 세대 전환
		
		// 다음 세대 적응도 평가
		nextGen_Fitness = determine_fitness(items, generation);

		/*/ 평가에 대한 결과 출력
		cout << gen_start + 2 << "번 세대 " << nextGen_Fitness.chromo_num + 1 <<
			"번째 유전자에서 최고 이익 : " << nextGen_Fitness.best_profit <<
			" 가방무게 : " << nextGen_Fitness.weight << " 변이 발생량 : " <<mutated_count <<
			" 선택변이 발생량 : " << selective_mutated_count << endl;
			*/
		gen_start += 1; // 세대 카운터
		//selective_mutated_count = 0; // 선택적 변이 횟수 카운터
		//mutated_count = 0; // 변이 횟수  카운터
		 
		// 각 라운드 통틀어 최고의 이익 갱신 (배낭의 무게를 넘지 않는 유전자만)
		if (nextGen_Fitness.weight < CAPACITY && bestChrom.best_profit < nextGen_Fitness.best_profit)
		{
			bestChrom.chromo_num = nextGen_Fitness.chromo_num;
			bestChrom.best_profit = nextGen_Fitness.best_profit;
			bestChrom.weight = nextGen_Fitness.weight;
			bestGen = gen_start + 1;
		}
		
	}
	// 결과발표
	cout << "최적의 유전자는 " << bestGen << "세대의 " << bestChrom.chromo_num + 1 << "번 유전자입니다.\n" << endl;
	cout << "최적 유전자의 profit은 " << bestChrom.best_profit << ", 무게는 " << bestChrom.weight << "입니다.\n";
	free(items);
}

void make_item()
{
	FILE* item_info = fopen("item_info.txt", "w");
	srand(time(NULL));
	int price = 0;
	int weight = 0;
	int k = 1;
	// 가방의 무게와 아이템 갯수를 고려하여 개별 아이템의 무게 범위를 구한다.
	// CAPACITY = 가방의무게, CL = 아이템 갯수
	// 아이템 무게는 1 ~ c_weight 에서 결정됨.
	// k는 정하고자 싶은 무게 total 무게범위 상수이다.
	// c_weight는 평균적으로 (1+k)*CAPACITY 가 나오게끔 구해진다.
	int c_weight = 4 * (1 + k)*(CAPACITY - 1) / CL;


	fprintf(item_info, "%d\n", CL);
	for (int i = 0; i < CL; i++)
	{
		price = 1 + (rand() % 100); // short 는 6만~
		weight = 1 + (rand() % c_weight); // weight 가 char이기 때문에 ~128까지
		fprintf(item_info, "%d %d\n", price, weight);
	}
	fclose(item_info);
}

void read_item(item *items)
{
	//FILE *fp = fopen("phasor.txt", "r");
	FILE *fp = fopen("item_info.txt", "r");
	int maxItem;
	int i_price;
	int i_weight;

	fscanf(fp, "%d", &maxItem);
	for (int i = 0; i < CL; i++)
	{
		fscanf(fp, "%d%d", &i_price, &i_weight);
		items[i].price = i_price;
		items[i].weight = i_weight;
	}
	fclose(fp);
}

// make first generation of chromosomes
void make_chromosomes()
{
	srand(time(NULL));
	int c_price = 0;
	int c_weight = 0;
	for (int i = 0; i < PS; i++)
	{
		for (int j = 0; j < CL; j++)
		{
			chromosomes[0][i][j] = rand() % 2; // 0, 1을 랜덤하게 부여하며 유전자 생성
		}
	}
	for (int k = 0; k < CL; k++)
	{
		fitness[k] = 1000;
	}
}

void print_chromosomes()
{
	for (int i = 0; i < PS; i++)
	{
		cout << i + 1 << ". ";
		for (int j = 0; j < CL; j++)
		{
			cout << chromosomes[0][i][j] << " ";
		}
		cout << endl;
	}
}

P determine_fitness(item *items, int generation)
{
	int length = 0;
	int c_weight;
	C chorm_cumulative[PS]; // 1. price 2. weight
	for (int i = 0; i < PS; i++) // 초기화
	{
		chorm_cumulative[i].first = 0;
		chorm_cumulative[i].second = 0;
	}
	int minimum_profit = CL * 100;
	P each_gen;


	for (int population = 0; population < PS; population++)
	{
		c_weight = 0;
		// 각 유전자마다의 누적 가치, 가격 계산
		for (length = 0; length < CL; length++)
		{
			if (chromosomes[generation][population][length]) // 만약 1이라면 items[j]를 배낭에 담는다.
			{
				chorm_cumulative[population].second += items[length].weight;
				chorm_cumulative[population].first += items[length].price;
			}
		}

		// 각 유전자의 최소 profit 갱신 
		if (minimum_profit > chorm_cumulative[population].first)
		{
			minimum_profit = chorm_cumulative[population].first;
		}
		// 각 유전자의 최고 profit 갱신해주면서 기록
		if (each_gen.best_profit < chorm_cumulative[population].first)
		{
			// 가장 뛰어난 유전자 정보 기억
			each_gen.best_profit = chorm_cumulative[population].first;
			each_gen.chromo_num = population;
			each_gen.weight = chorm_cumulative[population].second;
		}
	}

	if (each_gen.weight > CAPACITY)
	{
		// 선택적 적응 true로 하여 0이 더 많게끔 돌연변이를 일으킨다.
		each_gen.selective_mutate = true;
	}
	else
	{
		each_gen.selective_mutate = false;
	}

	// 적응도 평가
	for (int i = 0; i < PS; i++)
	{
		if (chorm_cumulative[i].second > CAPACITY) // 무게가 배낭 전체보다 커졌다면,
		{
			// 적응도 평가 가중치의 3/4으로
			fitness[i] = (3/4)*(chorm_cumulative[i].first - minimum_profit) + ((each_gen.best_profit - minimum_profit) / (FC - 1));
		}
		else
		{
			// 적응도 평가
			fitness[i] = (chorm_cumulative[i].first - minimum_profit) + ((each_gen.best_profit - minimum_profit) / (FC - 1));
		}
	}

	return each_gen;
}

void print_fitness(int generation)
{
	for (int population = 0; population < PS; population++)
	{
		cout << generation + 1 << "번째 세대 " << population + 1 << "번 적응도: " << fitness[population] << "\n";
	}
	cout << endl;
}

int rulet(int chromo)
{
	srand(time(NULL) + chromo);
	int j = 1, k = 0;
	double random = 0;
	double prob = 0; // i번째 유전자의 확률
	int sumFitness = 0;
	for (int i = 0; i < PS; i++)
	{
		sumFitness += fitness[i];
	}

	while (j)
	{
		random = (double)((double)(rand() % 100) / 100);
		for (k = 0; k < PS; k++)
		{
			prob = (double)((double)fitness[k] / (double)sumFitness);
			random -= prob;
			if (random <= 0)
			{
				j = 0;
				break;
			}
		}
	}

	return k;
}

void cross(int generation, int chromo1, int chromo2)
{
	// CP = 0.7 교차율
	int temp = 0;
	int crossPoint1 = CL * CP1;
	int crossPoint2 = CL * CP2;

	for (int i = crossPoint1; i < crossPoint2; i++)
	{
		temp = chromosomes[generation][chromo1][i];
		chromosomes[generation][chromo1][i] = chromosomes[generation][chromo2][i];
		chromosomes[generation][chromo2][i] = temp;
	}


}

void copy_chromosomes(int generation, int start, int chromo1, int chromo2)
{
	for (int j = 0; j < CL; j++)
	{
			chromosomes[(!generation)][start][j] = chromosomes[generation][chromo1][j];
			chromosomes[(!generation)][start][j+1] = chromosomes[generation][chromo2][j];
	}
}

void print_generation(int gen_start, int generation)
{
	cout << gen_start+2 << "세대의 유전자 정보.\n";
	for (int i = 0; i < PS; i++)
	{
		cout << gen_start + 2 << "-" << i + 1 << ". ";
		for (int j = 0; j < CL; j++)
		{

			cout << chromosomes[generation][i][j];
		}
		cout << "\n";
	}
}

void print_cross(int generation, int chromo1, int chromo2)
{
	int crossPoint = CL * CP1;
	cout << "croosPoint는 " << crossPoint << "입니다.\n";
	cout << "chromo1의 유전자 정보\n";
	for (int i = 0; i < CL; i++)
	{
		if (i == crossPoint)
			cout << " ";
		cout << chromosomes[generation][chromo1][i];
	}

	cout << "\nchromo2의 유전자 정보\n";
	for (int j = 0; j < CL; j++)
	{
		if (j == crossPoint)
			cout << " ";
		cout << chromosomes[generation][chromo2][j];
	}
	cout << "\n";
}

void mutation(int generation, int chromo1)
{
	//cout << generation + 1 << "번 세대 변이 시작\n";
	int random = 0;
	srand(time(NULL) + chromo1);
	if ((rand() % 100) + 1 <= (double)MP * 100) // 1% 확률로 변이
	{
		mutated_count += 1;
		for (int i = 0; i < MP * 50; i++) // 랜덤으로 10%의 단백질을 0으로 변이
		{
			random = rand() % 100;
			//랜덤으로 0~99까지의 단백질중 하나를 반전
			chromosomes[generation][chromo1][random] = !(chromosomes[generation][chromo1][random]);
			//cout << i + 1 << "유전자 " << random+1 << "번째 단백질\n";
		}
	}
}

void seletive_mutation(int generation, int chromo1)
{
	// 한 세대의 해가 capacity를 초과할 경우 의도적으로 아이템을 줄이는 변이를 시행.

	int random = 0;
	srand(time(NULL) + chromo1);

	if (rand() % 100 <= (double)SMP * 100) // 5% 확률로 변이
	{
		selective_mutated_count += 1;
		for (int i = 0; i < SMP * 2000; i++) // 랜덤으로 5%의 단백질을 0으로 변이
		{
			//랜덤으로 0~99까지의 단백질중 하나를
			random = rand() % 100;

			// 0으로 만든다.
			chromosomes[generation][chromo1][random] = 0;
			//cout << i + 1 << "유전자 " << random+1 << "번째 단백질\n";
		}
	}
	

}

int elitism_select(int start)
{
	srand(time(NULL) + start);
	int random = (rand() % 1000) * EP;

	return cpyFitness[random].first;
}

void insertion_sort(C *data, int n)
{
	int i, j;
	C remember;
	for (i = 1; i < n; i++)
	{
		remember = data[(j = i)];
		while (--j >= 0 && remember.second > data[j].second) {
			data[j + 1] = data[j];
			data[j] = remember;
		}
	}
}

void quickSort(C arr[], int left, int right)
{
	int i = left, j = right-1;
	int pivot = arr[(left + right) / 2].second;
	C temp;
	do
	{
		while (arr[i].second < pivot)
			i++;
		while (arr[j].second > pivot)
			j--;
		if (i <= j)
		{
			temp = arr[i];
			arr[i]= arr[j];
			arr[j] = temp;
			i++;
			j--;
		}
	} while (i <= j);
}

void copy_fitness()
{
	for (int i = 0; i < PS; i++)
	{
		cpyFitness[i].first = i;
		cpyFitness[i].second = fitness[i];
	}
}

void reNewChromosomes(int generation)
{
	int temp;
	srand(time(NULL));
	int random = rand() % CL;
	for (int i = 0; i < PS; i++)
	{
		for (int j = 0; i < CL; j++)
		{
			// random 자리와 j자리의 유전자를 교환
			temp = chromosomes[generation][i][j];
			chromosomes[generation][i][j] = chromosomes[generation][i][random];
			chromosomes[generation][i][random] = temp;
		}
	}
}

C make_randomNum()
{
	C randoms(101, 101);
	randoms.first = rand() % 1000;
	if (visited[randoms.first] == true)
	{
		while (visited[randoms.first])
		{
			randoms.first = rand() % 1000;
		}
	}
	visited[randoms.first] = true;


	randoms.second = rand() % 1000;
	if (visited[randoms.second] == true)
	{
		while (visited[randoms.second])
		{
			randoms.second = rand() % 1000;
		}
	}
	visited[randoms.second] = true;

	return randoms;
}

void init_visited()
{
	for (int i = 0; i < CL; i++)
	{
		visited[i] = false;
	}
}