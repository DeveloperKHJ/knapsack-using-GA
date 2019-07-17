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
#define CL 100 // chronomosome length ������ ����
#define PS 1000 // �� ������ �α�
#define GN 5000 // ������ ��(�ݺ� ��)
#define CAPACITY 3001// �ִ� ���빫��
#define CP1 0.3 // ����point1
#define CP2 0.7 // ���� point2
#define PC 0.9 // ������ �߻��ϴ� Ȯ�� 90%
#define MP 0.05 // ������ 1%
#define SMP 0.07 // ������ ������ (SMP*100)% �� �����ϰ� 0���� �����.
#define EP 0.04 // ����Ƽ�� ����, ���� (EP*100)%�� �������뿡 ���Խ�Ų��.
#define FC 3 // ������ ��� (k > 1) �������� ���þ� ����

// �������� ���� ����ü
struct item
{
	short weight;
	short price;
};

// �� ������, �� �����ڿ� ���� ������ ���� ����ü
struct P
{
	int chromo_num = 0;
	int best_profit = 0;
	int weight = 0;
	bool selective_mutate = false;
};

// chromosomes[����][�α�][������ ����]
// ���� : 0 -> ���� ����
//			  1 -> ���� ����
int chromosomes[2][PS][CL];


int selective_mutated_count = 0; // ���ú��� �߻�ī����
int mutated_count = 0; // ���� �߻�ī����
int fitness[PS]; // ������ �迭
bool visited[PS+1] = { false }; // ���� ������ �湮�迭, �����ϰ� �����Ҷ� �ߺ� ������
typedef pair<int, int> C; // 1. index, 2. fitness -> ������ ���Ľ� ���, in-place ������ ���� pair ���
C cpyFitness[PS]; // pair C�� �迭�� ���������ν� �������� in-place�ϰ� ���� ����

void make_item(); // item list�� �����
void read_item(item *items); // item list�� �ҷ��´�
void make_chromosomes(); // �ʱ� �����ڸ� �����Ѵ� �������� 0,1�� �ο��Ѵ�.
void print_chromosomes(); // �����ڸ� ����Ѵ�.
P determine_fitness(item *items, int generation); // �������� �����Ѵ�.
void print_fitness(int geneartion); // �������� ����Ѵ�.
int rulet(int chromo); // �귿 �Լ�
void cross(int generation, int chromo1, int chromo2); // chromo1�� chromo2�� ���� �Լ�
void print_cross(int generation, int chromo1, int chromo2); // ������ ����� ����Ѵ�.
void copy_chromosomes(int generation, int start, int chromo1, int chromo2); // �����ڸ� ��������� �����Ѵ�.
void print_generation(int gen_start, int generation); // ���� ��ü�� ����Ѵ�.
void mutation(int generation, int chromo1); // �����Լ�.
void seletive_mutation(int generation, int chromo1); // ������ ���� �Լ�
int elitism_select(int start); // ����Ƽ�� �Լ�
void copy_fitness(); // ������ �迭 ����
void insertion_sort(C *data, int n); // �������ķ� ����� ������ �迭�� ����
void init_visited(); // �湮�迭�� ��� false�� �ʱ�ȭ



int main()
{
	// item �迭
	int start = 0; // ������ �ݺ� Count
	int gen_start = 0; // ���� �ݺ� Count
	P nextGen_Fitness; // ������ �򰡿� ���
	P bestChrom; // ������ ������
	int bestGen = 0; // ������ �����ڰ� ���� ����
	int chromo1; // ���� ������ 1
	int chromo2; // ���� ������ 2 
	int generation = 0;
	item *items;
	items = (item *)malloc(sizeof(item)*CL);

	// ������ �ʱ�ȭ
	read_item(items);

	// ������ ����
	make_chromosomes();

	// �ʱ� ������ �˻�
	P init_gen = determine_fitness(items, generation);

	// ���� ����Ŭ ����
	while (gen_start < GN)
	{
		// ���� ���� ���� ����
		for (start = 0; start < PS - 1; start += 2)
		{
			// EP ������ŭ�� �����ڸ�
			// ����Ƽ������ ���� ���Ե� �����ڷ� ����
			if (start < EP*PS)
			{
				copy_fitness(); // �������� ���� �ѵ�
				insertion_sort(cpyFitness, PS); // ������ ����
				chromo1 = elitism_select(start + 100); // �������� ���� ����
				chromo2 = elitism_select(start + 50000); // ���ڰ��� random seed
			}
			else // ����Ƽ�� ������ ������ ������ �귿���� ����
			{
				chromo1 = rulet(start + gen_start+100); // ������ �迭�� ������� �� �귿�Լ�
				chromo2 = rulet(start); // ���ڰ��� random seed
			}

			// ����
			srand(time(NULL) + gen_start + start); // random seed
			if (rand() % 100 <= PC * 100)
			{
				cross(generation, chromo1, chromo2); // ���õ� chromo1 �� chromo2�� ����
			}
			
			// ���� ���� ���� ��, ���� ���뿡 ���� ����
			// ������ ���Ը� �Ѵ� �����ڿ� ���ؼ��� ������ ���̸� ����
			if (nextGen_Fitness.selective_mutate)
			{
				seletive_mutation(generation, start); // ������ ���̴� ������ �ƴ϶� 0���� �����.
				seletive_mutation(generation, start + 1); // start�� 2���� �ڹǷ� 2����
			}
			else
			{
				// ������ ���Ը� ���� �ʴ� ��� ���� ����
				mutation(generation, start);
				mutation(generation, start + 1);
			}

			// �������뿡 ����
			copy_chromosomes(generation, start, chromo1, chromo2);
		}
		init_visited(); // �湮 �迭 false�� �ʱ�ȭ
		generation = !generation; // ���� ���� ��ȯ
		
		// ���� ���� ������ ��
		nextGen_Fitness = determine_fitness(items, generation);

		/*/ �򰡿� ���� ��� ���
		cout << gen_start + 2 << "�� ���� " << nextGen_Fitness.chromo_num + 1 <<
			"��° �����ڿ��� �ְ� ���� : " << nextGen_Fitness.best_profit <<
			" ���湫�� : " << nextGen_Fitness.weight << " ���� �߻��� : " <<mutated_count <<
			" ���ú��� �߻��� : " << selective_mutated_count << endl;
			*/
		gen_start += 1; // ���� ī����
		//selective_mutated_count = 0; // ������ ���� Ƚ�� ī����
		//mutated_count = 0; // ���� Ƚ��  ī����
		 
		// �� ���� ��Ʋ�� �ְ��� ���� ���� (�賶�� ���Ը� ���� �ʴ� �����ڸ�)
		if (nextGen_Fitness.weight < CAPACITY && bestChrom.best_profit < nextGen_Fitness.best_profit)
		{
			bestChrom.chromo_num = nextGen_Fitness.chromo_num;
			bestChrom.best_profit = nextGen_Fitness.best_profit;
			bestChrom.weight = nextGen_Fitness.weight;
			bestGen = gen_start + 1;
		}
		
	}
	// �����ǥ
	cout << "������ �����ڴ� " << bestGen << "������ " << bestChrom.chromo_num + 1 << "�� �������Դϴ�.\n" << endl;
	cout << "���� �������� profit�� " << bestChrom.best_profit << ", ���Դ� " << bestChrom.weight << "�Դϴ�.\n";
	free(items);
}

void make_item()
{
	FILE* item_info = fopen("item_info.txt", "w");
	srand(time(NULL));
	int price = 0;
	int weight = 0;
	int k = 1;
	// ������ ���Կ� ������ ������ ����Ͽ� ���� �������� ���� ������ ���Ѵ�.
	// CAPACITY = �����ǹ���, CL = ������ ����
	// ������ ���Դ� 1 ~ c_weight ���� ������.
	// k�� ���ϰ��� ���� ���� total ���Թ��� ����̴�.
	// c_weight�� ��������� (1+k)*CAPACITY �� �����Բ� ��������.
	int c_weight = 4 * (1 + k)*(CAPACITY - 1) / CL;


	fprintf(item_info, "%d\n", CL);
	for (int i = 0; i < CL; i++)
	{
		price = 1 + (rand() % 100); // short �� 6��~
		weight = 1 + (rand() % c_weight); // weight �� char�̱� ������ ~128����
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
			chromosomes[0][i][j] = rand() % 2; // 0, 1�� �����ϰ� �ο��ϸ� ������ ����
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
	for (int i = 0; i < PS; i++) // �ʱ�ȭ
	{
		chorm_cumulative[i].first = 0;
		chorm_cumulative[i].second = 0;
	}
	int minimum_profit = CL * 100;
	P each_gen;


	for (int population = 0; population < PS; population++)
	{
		c_weight = 0;
		// �� �����ڸ����� ���� ��ġ, ���� ���
		for (length = 0; length < CL; length++)
		{
			if (chromosomes[generation][population][length]) // ���� 1�̶�� items[j]�� �賶�� ��´�.
			{
				chorm_cumulative[population].second += items[length].weight;
				chorm_cumulative[population].first += items[length].price;
			}
		}

		// �� �������� �ּ� profit ���� 
		if (minimum_profit > chorm_cumulative[population].first)
		{
			minimum_profit = chorm_cumulative[population].first;
		}
		// �� �������� �ְ� profit �������ָ鼭 ���
		if (each_gen.best_profit < chorm_cumulative[population].first)
		{
			// ���� �پ ������ ���� ���
			each_gen.best_profit = chorm_cumulative[population].first;
			each_gen.chromo_num = population;
			each_gen.weight = chorm_cumulative[population].second;
		}
	}

	if (each_gen.weight > CAPACITY)
	{
		// ������ ���� true�� �Ͽ� 0�� �� ���Բ� �������̸� ����Ų��.
		each_gen.selective_mutate = true;
	}
	else
	{
		each_gen.selective_mutate = false;
	}

	// ������ ��
	for (int i = 0; i < PS; i++)
	{
		if (chorm_cumulative[i].second > CAPACITY) // ���԰� �賶 ��ü���� Ŀ���ٸ�,
		{
			// ������ �� ����ġ�� 3/4����
			fitness[i] = (3/4)*(chorm_cumulative[i].first - minimum_profit) + ((each_gen.best_profit - minimum_profit) / (FC - 1));
		}
		else
		{
			// ������ ��
			fitness[i] = (chorm_cumulative[i].first - minimum_profit) + ((each_gen.best_profit - minimum_profit) / (FC - 1));
		}
	}

	return each_gen;
}

void print_fitness(int generation)
{
	for (int population = 0; population < PS; population++)
	{
		cout << generation + 1 << "��° ���� " << population + 1 << "�� ������: " << fitness[population] << "\n";
	}
	cout << endl;
}

int rulet(int chromo)
{
	srand(time(NULL) + chromo);
	int j = 1, k = 0;
	double random = 0;
	double prob = 0; // i��° �������� Ȯ��
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
	// CP = 0.7 ������
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
	cout << gen_start+2 << "������ ������ ����.\n";
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
	cout << "croosPoint�� " << crossPoint << "�Դϴ�.\n";
	cout << "chromo1�� ������ ����\n";
	for (int i = 0; i < CL; i++)
	{
		if (i == crossPoint)
			cout << " ";
		cout << chromosomes[generation][chromo1][i];
	}

	cout << "\nchromo2�� ������ ����\n";
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
	//cout << generation + 1 << "�� ���� ���� ����\n";
	int random = 0;
	srand(time(NULL) + chromo1);
	if ((rand() % 100) + 1 <= (double)MP * 100) // 1% Ȯ���� ����
	{
		mutated_count += 1;
		for (int i = 0; i < MP * 50; i++) // �������� 10%�� �ܹ����� 0���� ����
		{
			random = rand() % 100;
			//�������� 0~99������ �ܹ����� �ϳ��� ����
			chromosomes[generation][chromo1][random] = !(chromosomes[generation][chromo1][random]);
			//cout << i + 1 << "������ " << random+1 << "��° �ܹ���\n";
		}
	}
}

void seletive_mutation(int generation, int chromo1)
{
	// �� ������ �ذ� capacity�� �ʰ��� ��� �ǵ������� �������� ���̴� ���̸� ����.

	int random = 0;
	srand(time(NULL) + chromo1);

	if (rand() % 100 <= (double)SMP * 100) // 5% Ȯ���� ����
	{
		selective_mutated_count += 1;
		for (int i = 0; i < SMP * 2000; i++) // �������� 5%�� �ܹ����� 0���� ����
		{
			//�������� 0~99������ �ܹ����� �ϳ���
			random = rand() % 100;

			// 0���� �����.
			chromosomes[generation][chromo1][random] = 0;
			//cout << i + 1 << "������ " << random+1 << "��° �ܹ���\n";
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
			// random �ڸ��� j�ڸ��� �����ڸ� ��ȯ
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