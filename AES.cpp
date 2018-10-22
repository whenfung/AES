#include "AES_tables.h"
#include <string.h>

unsigned char w[44][4];           //10个轮密钥

unsigned char SBox_Chg(char sin,int bits) {  //对字符sin进行字节替换 
//S盒变换
	int i,j,k=0;
	char temp,out[8];
	for(j=7; j>=0; j--) { //获得行号列号
		temp=sin&(1<<j);
		if(0 == temp) {
			out[k]=0;
			k++;
		} else {
			out[k]=1;
			k++;
		}
	}
	i=out[0]*8+out[1]*4+out[2]*2+out[3];  //行号
	j=out[4]*8+out[5]*4+out[6]*2+out[7];  //列号
	return (S_box[i][j]);   //返回得到的S盒的值
}

unsigned char SBox_1_Chg(char sin,int bits) {  //就返回值和S盒替换不一样 
//S盒逆变换
	int i,j,k=0;
	char temp,out[8];
	for(j=7; j>=0; j--) {
		temp=sin&(1<<j);
		if(0 == temp) {
			out[k]=0;
			k++;
		} else {
			out[k]=1;
			k++;
		}
	}
	i=out[0]*8+out[1]*4+out[2]*2+out[3];//行号
	j=out[4]*8+out[5]*4+out[6]*2+out[7];//列号
	return (SBox_1[i][j]);      

}

void Getkey(unsigned char key[16]) {  
	//10个轮密钥生成
	static unsigned char temp[4],t[4];
	int i,j,k;
	for(i=0; i<4; i++) {       
		memcpy(w[i],&key[i*4],4);     //这里的轮密钥按行存储，每行32bits，0到4行是原始密钥 
	}
	for(i=4; i<44; i++) {               
		memcpy(temp,w[i-1],4);          //w[j-1] 
		if(i%4==0) {                    //if(j-1)%4==0,进行变换 
			memcpy(t,temp,1);           //3个memcpy其实是进行循环左移一位 
			memcpy(temp,&temp[1],3);
			memcpy(&temp[3],t,1);
			for(j = 0,k=(i/4)*4; j<4; j++,k++) {          //进行S盒变换后和RC进行异或 
				temp[j]=SBox_Chg(temp[j],8)^RC[k-4];      //每次循环和RC数组异或一个字节，循环四次，生成g(w[j-1]) 
			}
			for(j=0; j<4; j++) {                          //w[j] = g(w[j-1]) ^ w[i-4] 
				w[i][j]=w[i-4][j]^temp[j];
			}
		} else
			for(j=0; j<4; j++) {                         //w[j] = w[j-4] ^ w[j-1]
				w[i][j]=w[i-4][j]^temp[j];
			}
	}
}

void rows_mov(unsigned char *sin) {  //行移位运算
	static unsigned char temp[4];
	memcpy(temp,&sin[4],1);     //第一行不变，第二行循环左移一位
	memcpy(&sin[4],&sin[5],3);
	memcpy(&sin[7],temp,1);

	memcpy(temp,&sin[8],2);      //第三行循环左移两位
	memcpy(&sin[8],&sin[10],2);
	memcpy(&sin[10],temp,2);

	memcpy(temp,&sin[12],3);     //第四行循环左移三位
	memcpy(&sin[12],&sin[15],1);
	memcpy(&sin[13],temp,3);
}

void rows_mov_1(unsigned char *sin) {  //行移位运算逆运算
	static unsigned char temp[4];
	memcpy(temp,&sin[7],1);
	memcpy(&sin[5],&sin[4],3);
	memcpy(&sin[4],temp,1);       //循环右移1位，实质为内存内容相互交换

	memcpy(temp,&sin[10],2);
	memcpy(&sin[10],&sin[8],2);
	memcpy(&sin[8],temp,2);      //循环右移2位，实质为内存内容相互交换

	memcpy(temp,&sin[12],1);
	memcpy(&sin[12],&sin[13],3);
	memcpy(&sin[15],temp,1);      //循环右移3位，实质为内存内容相互交换
}

unsigned char GF2mul(unsigned char a, unsigned char b) { 
//有限域GF(2^8)上的乘法
	unsigned char bw[4];
	unsigned char res=0;
	int i;
	bw[0] = b;
	for(i=1; i<4; i++) {       //循环三次，得到乘2、4、8后的值
		bw[i] = bw[i-1]<<1;    //原数值乘2
		if(bw[i-1]&0x80) {     //判断原数值是否小于0x80，掩码原理 
			bw[i]^=0x1b;       //最左为1，减去0x1b doublesand 
		}
	}
	for(i=0; i<4; i++) {
		if((a>>i)&0x01) {    //将参数a的值表示为1、2、4、8的线性组合
			res ^= bw[i];
		}
	}
	return res;
}

void columnsmix(unsigned char *sin) {
//列混合变换
	unsigned char t[4][4]= {0};
	int i,j,k;
	for(i=0; i<4; i++)
		for(k=0; k<4; k++)
			for(j=0; j<4; j++)
				t[i][k]^=GF2mul(L_mix[i][j],sin[j*4+k]);  //矩阵运算：c[i][k] += a[i][j] * b[j][k]
	memcpy(sin,t,16);
}


void columnsmix_1(unsigned char *sin) {
//逆列混合变换
	unsigned char t[4][4]= {0};
	int i,j,k;
	for(i=0; i<4; i++)
		for(k=0; k<4; k++)
			for(j=0; j<4; j++)
				t[i][k]^=GF2mul(L_mix_1[i][j],sin[j*4+k]);
	memcpy(sin,t,16);
}

void addroundkey_start(unsigned char *sin,unsigned char (*p)[4]) {  //sin[j][i] = sin[i][j]^p[i][j]
//最初一次轮密钥加操作，按顺序异或 
	unsigned char out[16];
	int i,j,k=0;
	for(i=0; i<4; i++)
		for(j=0; j<4; j++) {
			out[j*4+i]=sin[k]^p[i][j];  
			k++;
		}
	memcpy(sin,out,16);
}


void addroundkey(unsigned char *sin,unsigned char (*p)[4]) {
//轮密钥加操作，相应的为进行异或，由于第一次轮密钥之后sin逻辑上是按列排列的，即 sin[j][i] = sin[j][i]^p[i][j] 
	unsigned char out[16];
	int i,j,k=0;
	for(i=0; i<4; i++)
		for(j=0; j<4; j++) {
			out[j*4+i]=sin[j*4+i]^p[i][j];
			k++;
		}
	memcpy(sin,out,16);
}

void AES_Cry(unsigned char *sin) {   //sin是信息的意思 
	int i,j,k=0;

	addroundkey_start(sin,&w[0]);   //第一轮密钥加操作 
	for(i=1; i<10; i++) {            //S盒置换 
		for(j=0; j<16; j++)
			sin[j]=SBox_Chg(sin[j],8);
		rows_mov(sin);                   //行移位 
		columnsmix(sin);                 //列混淆 
		addroundkey(sin,&w[i*4]);        //轮密钥加 

	}
	for(j=0; j<16; j++)                //第十六轮加密  
		sin[j]=SBox_Chg(sin[j],8);      //S盒置换 
	rows_mov(sin);                      //行移位 
	addroundkey(sin,&w[i*4]);           //轮密钥加 
}

void AES_Dec(unsigned char *sin) {
	int i,j,k=0;                     

	addroundkey(sin,&w[40]);           //轮密钥加 
	for(i=9; i>0; i--) {              
		rows_mov_1(sin);                    //逆向行移位  
		for(j=0; j<16; j++)
			sin[j]=SBox_1_Chg(sin[j],8);    //逆向字节替换 
		addroundkey(sin,&w[i*4]);           //轮密钥加 
		columnsmix_1(sin);                  //逆向列混淆 
	}
	rows_mov_1(sin);                        //逆向行移位 
	for(j=0; j<16; j++)                    //逆向字节替换 
		sin[j]=SBox_1_Chg(sin[j],8);        
	addroundkey(sin,&w[i*4]);               //轮密钥加 
}

int main() {
	int i,j;

	unsigned char messages[17]= {0}; //不足16个字符自动补0 
	unsigned char Mykey[17]= {0};    
	unsigned char Yourkey[17]= {0};  
	printf("输入明文:(128bit)\n");
	gets((char*)messages);  //输入messages的时候是char型，不影响后面的加密 

	printf("输入密钥:(128bit)\n");
	gets((char*)Mykey);

	Getkey(Mykey);
	AES_Cry(messages);
	printf("\n\n加密后：");
	for(i=0; i<4; i++)
		for(j=0; j<4; j++)
			printf("%x ",messages[j*4+i]);              
	printf("\n输入密钥:(128bit)\n");
	gets((char*)Yourkey);
	Getkey(Yourkey);
	AES_Dec(messages);
	printf("\n解密后：");
	for(i=0; i<4; i++)   //解密后的字符串不是顺序的，实际是按列存储的4*4字符矩阵 
		for(j=0; j<4; j++)
			printf("%c ",messages[j*4+i]);
	printf("\n");
}
