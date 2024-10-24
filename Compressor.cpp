#include <stdio.h>
#include <vector>
#include <map>
#include <iostream> 
#include <sstream>
#include <fstream>
#include <string>
#include <numeric>
#include <time.h>
#include <chrono>
#include <thread>
#include <cstring>
#include <cstdlib>
#include <unistd.h>

// TODO usar 4096
#define K 16
#define SIZE_DICIONARIO 65536
#define SIZE 2048
#define arqEntrada "SRR22399485_1.fastq"
#define arqSaida "SRR22399485_1.lzw"
#define arqDecodificado "decodificado.fastq"


using namespace std;
std::map<std::vector<int>, int> dicionario;
std::map<int, std::vector<int> > dicionarioInverso;
unsigned char* buffer = new unsigned char[SIZE];
int pos = 0; // posicao em bits
FILE* output;
//const char *name = "output.lzw";
const char* nomeArquivoGene = "gene.txt";
const char* nomeArquivoGeral = "geral.txt";
int tamanhoCod;
int codigo = 128;
std::vector<int> texto;
std::vector<int> sequencia;
std::vector<int> sequenciaDeSaida;

void addBits(int value, int qBits)
{
    int byte = pos / 8; // posicao do proximo byte com bits disponiveis
    int shift = pos % 8; // quant de bits utilizados
    int disp = 8 - shift; // quant de bits disponiveis
    unsigned int mask = 0;
    int falta = qBits - disp; // bits que faltam para o proximo byte
    int falta2 = 0; // bits que faltam para o segundo proximo byte
    int i;
    if (falta > 8) {
        falta2 = falta - 8;
        falta = 8;
    }

    //printf("value=%3d::byte=%d, shift=%d, disp=%d, mask=%d, falta=%d, falta2=%d\n", value, byte, shift, disp, mask, falta, falta2);

    mask = 0;
    for (i = 0; i < disp; i++)
        mask = mask | (1 << i);
    buffer[byte] = (buffer[byte] & (~mask & 0xFF)) | ((value >> (falta + falta2)) & mask);
    //printf("mask1=0x%02X (~0x%02X)\n", mask, (~mask & 0xFF));

    if (falta > 0) {
        mask = 0;
        for (i = 0; i < (8 - falta); i++)
            mask = mask | (1 << i);
        if (falta2 == 0)
            buffer[byte + 1] = (buffer[byte + 1] & mask) | ((value << (8 - falta)) & (~mask & 0xFF));
        else
            buffer[byte + 1] = ((value >> falta2) & 0xFF);
        //printf("mask2=0x%02X (~0x%02X)\n", mask, (~mask & 0xFF));
    }

    if (falta2 > 0) {
        mask = 0;
        for (i = 0; i < (8 - falta2); i++)
            mask = mask | (1 << i);
        buffer[byte + 2] = (buffer[byte + 2] & mask) | ((value << (8 - falta2)) & (~mask & 0xFF));
        //printf("mask3=0x%02X (~0x%02X)\n", mask, (~mask & 0xFF));
    }

    pos += qBits;

    // Adicionado para salvar os bytes as poucos.
    byte = pos / 8; // posicao do proximo byte com bits disponiveis
    if (byte >= SIZE - 3) {
        fwrite(buffer, 1, byte, output);
        //printf("salvei %d bytes\n", byte);
        pos -= byte * 8;
        shift = pos % 8; // quant de bits utilizados
        for (i = 0; i < (8 - shift); i++)
            mask = mask | (1 << i);
        buffer[0] = buffer[byte] & ~mask;
    }
}

int getBits(int qBits)
{
    int value = 0;
    int byte = pos / 8; // posicao do proximo byte com bits disponiveis
    int shift = pos % 8; // quant de bits utilizados
    int disp = 8 - shift; // quant de bits disponiveis
    unsigned int mask = 0;
    int falta = qBits - disp; // bits que faltam para o proximo byte
    int falta2 = 0; // bits que faltan para o segundo proximo byte
    int i;
    if (falta > 8) {
        falta2 = falta - 8;
        falta = 8;
    }

    //printf("value=%3d::byte=%d, shift=%d, disp=%d, mask=%d, falta=%d, falta2=%d\n", value, byte, shift, disp, mask, falta, falta2);

    mask = 0;
    for (i = 0; i < disp; i++)
        mask = mask | (1 << i);
    value = (buffer[byte] & mask) << (falta + falta2);
    //printf("mask1=0x%02X (~0x%02X)\n", mask, (~mask & 0xFF));

    if (falta > 0) {
        if (falta2 == 0)
            value = value | ((buffer[byte + 1]) >> (8 - falta));
        else
            value = value | (buffer[byte + 1] << falta2);
        //printf("mask2=0x%02X (~0x%02X)\n", mask, (~mask & 0xFF));
    }

    if (falta2 > 0) {
        //mask = 0;
        //for (i = 0 ; i < (8 - falta2) ; i++)
        //    mask = mask | (1 << i);
        value = value | (buffer[byte + 2] >> (8 - falta2));
        //printf("mask3=0x%02X (~0x%02X)\n", mask, (~mask & 0xFF));
    }
    pos += qBits;
    return value;
}

void salvaFim() {
    int byte = pos / 8; // posicao do proximo byte com bits disponiveis
    int shift = pos % 8;
    if (shift > 0)
        byte++;
    if (byte > 0) {
        fwrite(buffer, 1, byte, output);
        //printf("salvei fim %d bytes\n", byte);
    }
    fclose(output);
}

void lerArquivo(const char* name) {
    texto.clear();
    unsigned char byte;
    FILE* input = fopen(name, "rb");
    if (input == NULL) {
        printf("Erro abrindo o arquivo %s\n", name);
        return;
    }

    while (!feof(input)) {
        size_t tam = fread(&byte, 1, 1, input);
        if (feof(input))
            break;
        //printf("%02X ", byte);
        texto.push_back(byte);
    }
    fclose(input);
}

bool verificaPresenca(std::vector<int> palavraAtual) {

    if (dicionario.count(palavraAtual) == 0) return false;
    return true;
}

bool verificaPresencaUndoer(int cod) {

    if (dicionarioInverso[cod].empty()) return false;
    else return true;
}

int execCommand(const char* command)
{
    int result = system(command);
    if (result == 0)
    {
        return -1;
    }
    else
    {
        return 1;
    }
}

void inicializaDicionario()
{
    for (int i = 0; i < 128; i++) {
        dicionario[{i}] = i;
    }
}

void inicializaDicionarioInverso()
{
    for (int i = 0; i < 128; i++) {
        dicionarioInverso[i] = { i };
    }
}

void splitterLzw()
{
    bool usandoGene = false;
    int countLine = 3;
    ofstream arqGene, arqGeral;
    arqGene.open(nomeArquivoGene, ios::trunc);
    arqGeral.open(nomeArquivoGeral, ios::trunc);
    if (arqGeral.is_open() && arqGene.is_open())
    {
        for (char i : texto)
        {
            if (usandoGene)
            {
                arqGene << i;
                if (i == '\n')
                {
                    usandoGene = false;
                    countLine++;
                }
            }
            else
            {
                arqGeral << i;
                if (i == '\n')
                {
                    countLine++;
                    if (countLine == 4)
                    {
                        usandoGene = true;
                        countLine = 0;
                    }
                }
            }

        }
    }
    arqGene.close();
    arqGeral.close();
}

void zipper()
{
    ifstream arqGene(nomeArquivoGene), arqGeral(nomeArquivoGeral);
    ofstream arqDec;
    arqDec.open(arqDecodificado, ios::trunc);
    int contLinhas = 3;
    string linha;
    while (getline(arqGeral, linha))
    {
        arqDec << linha;
        arqDec << '\n';
		contLinhas++;
        if (contLinhas == 4)
        {
        	getline(arqGene, linha);
        	arqDec << linha;
        	arqDec << '\n';
        	contLinhas=1;
        }
    }
    arqDec.close();
    arqGene.close();
    arqGeral.close();
}

void lzw(const char* arqToCod)
{
    lerArquivo(arqToCod);
    int qBits = K;
    
    //Percorre o texto
    std::vector<int> palavraAtual = { texto[0] };
    for (int i = 1; i < texto.size(); i = i) {

        while (verificaPresenca(palavraAtual) && i < texto.size()) {
            palavraAtual.push_back(texto[i]);
            i++;
        }

        if (!verificaPresenca(palavraAtual)) {
            if (dicionario.size() < SIZE_DICIONARIO) {

                dicionario[palavraAtual] = codigo++;
            }
            palavraAtual.pop_back();
            i--;
        }
        addBits(dicionario[palavraAtual], qBits);

        sequencia.push_back(dicionario[palavraAtual]);
        palavraAtual = { texto[i++] };

        if (i == texto.size()) {
            addBits(dicionario[palavraAtual], qBits);
            sequencia.push_back(dicionario[palavraAtual]);
        }

    }
    addBits(0, qBits);
    sequencia.clear();
    
}

void lzwUndoer(const char* arqToWrite)
{
    output = fopen(arqToWrite, "wb");
    if (output == NULL) {
        printf("Erro abrindo o arquivo %s\n", arqToWrite);
        return;
    }

    codigo = 128;
    int codigoAtual = texto[0];

    std::vector<int> sequenciaAnterior = { texto[0] };
    std::vector<int> sequenciaAtual = { texto[0] };

    fwrite(&dicionarioInverso[codigoAtual][0], 1, 1, output);
    sequenciaDeSaida.push_back(dicionarioInverso[codigoAtual][0]);

    for (int i = 1; i < texto.size(); i++) {

        codigoAtual = texto[i];
        sequenciaAnterior = sequenciaAtual;

        if (verificaPresencaUndoer(codigoAtual)) {
            // Código existe no dicionário
            sequenciaAtual = dicionarioInverso[codigoAtual];

            for (int j = 0; j < sequenciaAtual.size(); j++) {

                fwrite(&sequenciaAtual[j], 1, 1, output);
                sequenciaDeSaida.push_back(sequenciaAtual[j]);
            }

            if (codigo < SIZE_DICIONARIO) {
                for (int j = 0; j < sequenciaAnterior.size(); j++) {
                    dicionarioInverso[codigo].push_back(sequenciaAnterior[j]);
                }
                dicionarioInverso[codigo].push_back(sequenciaAtual[0]);
                codigo++;
            }
        }
        else {
            // Código não existe no dicionário

            for (int j = 0; j < sequenciaAnterior.size(); j++) {
                if (codigo < SIZE_DICIONARIO) {

                    dicionarioInverso[codigo].push_back(sequenciaAnterior[j]);
                }
                fwrite(&sequenciaAnterior[j], 1, 1, output);
                sequenciaDeSaida.push_back(sequenciaAnterior[j]);
            }

            if (codigo < SIZE_DICIONARIO) {

                dicionarioInverso[codigo].push_back(sequenciaAnterior[0]);
                codigo++;
            }
            fwrite(&sequenciaAnterior[0], 1, 1, output);
            sequenciaDeSaida.push_back(sequenciaAnterior[0]);
            sequenciaAtual.push_back(sequenciaAnterior[0]);

        }
    }
    fclose(output);
}

void splitterLzwUndoer()
{
    output = fopen(arqSaida, "rb");
    if (output == NULL) {
        printf("Erro abrindo o arquivo %s\n", arqSaida);
        return;
    }

    fseek(output, 0, SEEK_END);
    int tam = ftell(output);
    fseek(output, 0, SEEK_SET);
    if (buffer)
        delete[]buffer;
    buffer = new unsigned char[tam];
    
    size_t lange = fread(buffer, 1, tam, output);
    int cont = 0;
    fclose(output);
    for (int i = 0; i < tam / 2; i++)
    {
        int caracter = getBits(K);
        if (caracter == 0)
        {
            lzwUndoer((cont == 0) ? nomeArquivoGeral : nomeArquivoGene);
            texto.clear();
            dicionarioInverso.clear();
            sequenciaDeSaida.clear();
            inicializaDicionarioInverso();
            if (cont == 0) { cont = 1; continue; }
            else break;
        }
        texto.push_back(caracter);
    }

}

void callLzw()
{
	char command[100];
	sprintf(command, "touch %s %s %s", nomeArquivoGene, nomeArquivoGeral, arqSaida);
	execCommand(command);
	sprintf(command, "truncate -s 0 %s", arqSaida);
	execCommand(command);
	this_thread::sleep_for(std::chrono::milliseconds(100));
	output = fopen(arqSaida, "ab");
    if (output == NULL) {
        printf("Erro abrindo o arquivo %s\n", arqSaida);
        return;
    }
    inicializaDicionario();
	lerArquivo(arqEntrada);
	splitterLzw();
    /*lzw(nomeArquivoGeral);
    dicionario.clear();
    inicializaDicionario();
    codigo=128;
    lzw(nomeArquivoGene);
    sprintf(command, "rm %s %s", nomeArquivoGene, nomeArquivoGeral);
    execCommand(command);
    salvaFim();   */
}

void callLzwUndoer()
{
	char command[100];
	sprintf(command, "touch %s %s %s", nomeArquivoGene, nomeArquivoGeral, arqDecodificado);
	execCommand(command);
	sprintf(command, "truncate -s 0  %s %s %s", arqDecodificado, nomeArquivoGene, nomeArquivoGeral);
	execCommand(command);
	this_thread::sleep_for(std::chrono::milliseconds(100));
	inicializaDicionarioInverso();
	splitterLzwUndoer();
	zipper();
	sprintf(command, "rm %s %s", nomeArquivoGene, nomeArquivoGeral);
	execCommand(command);
}

int main()
{		clock_t begin = clock();
		callLzw();
	//	callLzwUndoer();
	  clock_t end = clock();
    cout << "Tempo LZW = " << (double)(end - begin) / CLOCKS_PER_SEC << endl;

}
