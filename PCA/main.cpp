#include"gdal_priv.h"
#include"cpl_conv.h" // for CPLMalloc()
#include<iostream>
#include<fstream>
#include<string>
#include"Eigen/Dense"
using namespace std;
using namespace Eigen;

//eigenʵ�����ɷַ���
void featurenormalize(MatrixXd &X)
{
	//����ÿһά�Ⱦ�ֵ
	MatrixXd meanval = X.colwise().mean();
	RowVectorXd meanvecRow = meanval;
	//������ֵ��Ϊ0
	X.rowwise() -= meanvecRow;
}
void computeCov(MatrixXd &X, MatrixXd &C)
{
	//����Э�������C = XTX / n-1;
	C = X.adjoint() * X;
	C = C.array() / (X.rows() - 1);
}
void computeEig(MatrixXd &C, MatrixXd &vec, MatrixXd &val)
{
	//��������ֵ������������ʹ��selfadjont���ն��������㷨ȥ���㣬�����ò�����vec��val������������
	SelfAdjointEigenSolver<MatrixXd> eig(C);

	vec = eig.eigenvectors();
	val = eig.eigenvalues();

}
int computeDim(MatrixXd &val)
{

	double sum = 0;
	for (int i = val.rows() - 1; i >= 0; --i)
	{
		sum = val(i, 0)/ val.sum();
	

		cout << "��"<<7-i<<"���ɷ�"<<"������:" << sum << endl;
			
	}
	/*int dim;
	double sum = 0;
	for (int i = val.rows() - 1; i >= 0; --i)
	{
	sum += val(i, 0);
	dim = i;

	if (sum / val.sum() >= 0.95)
	break;
	}
	return val.rows() - dim;*/
	return 7;
}

void writePcaImg(const char* path, int width, int height, double *pBuff, double *adfGeo, const char *prj, int bandNum, int imageSize, int pcaInd)
{
	GDALDriver *pDriver = GetGDALDriverManager()->GetDriverByName("GTiff"); //ͼ������
	char** ppszOptions = NULL;
	int depth = 8;//ͼ��λ��
	int dim = 1;//ÿ��ͼ�񲨶��������ｫÿ�����ɷִ洢��һ��������ͼ��
	GDALDataset* dst = pDriver->Create(path, width, height, dim, GDT_Float64, ppszOptions);//����ͼ��
	if (dst == nullptr)
		printf("Can't Write Image!");
	dst->SetGeoTransform(adfGeo);//��������
	dst->SetProjection(prj);//����ͶӰ
	dst->RasterIO(GF_Write, 0, 0, width, height, &pBuff[(bandNum - pcaInd)*imageSize], width, height,
		GDT_Float64, dim, nullptr, dim*depth, width*dim*depth, depth);//д��ͼ��
	GDALClose(dst);
}
int main(int argc, char *argv[])
{	//��ȡӰ��
	char* pszFilename = "D:/gdalData/pca/before.img";
	char *outPath = "D:/pca_temp/pca";
	/*char* pszFilename = argv[1];
	char* outPath = argv[2];*/
	GDALDataset  *poDataset;
	GDALAllRegister();
	poDataset = (GDALDataset *)GDALOpen(pszFilename, GA_ReadOnly);
	if (poDataset == NULL)
	{
		printf_s("read failed!\n");
	}
	else
	{
		printf_s("read successful!\n");
	}
	double adfGeoTransform[6];
	if (poDataset->GetGeoTransform(adfGeoTransform) == CE_Failure)//��ȡ������Ϣ
	{
		printf("��ȡ����ʧ��");
	}
	const char *prj = poDataset->GetProjectionRef();//��ȡͶӰ��Ϣ


	int iWidth = poDataset->GetRasterXSize();//ͼ����
	int iHeight = poDataset->GetRasterYSize();//ͼ��߶�
	int iBandCount = poDataset->GetRasterCount();//������
	int iImageSize = iWidth * iHeight;//ͼ����Ԫ��

	double *pBuff1 = new double[iImageSize*iBandCount];//���ٿռ�洢ԭʼͼ��

	poDataset->RasterIO(GF_Read, 0, 0, iWidth, iHeight, pBuff1,
		iWidth, iHeight, GDT_Float64, iBandCount, 0, 0, 0, 0);//��ȡԭʼͼ��

	MatrixXd staMat = Map<MatrixXd>(pBuff1, iImageSize, iBandCount);//��ͼ�����eigen����


	MatrixXd X(iImageSize, iBandCount), C(iBandCount, iBandCount);//�����δ洢��X���󣬹���Э�������C
	MatrixXd vec, val;//������������������ֵ����vec��val


	X = MatrixXd(staMat);

	//���ֵ��
	featurenormalize(X);

	//����Э����
	computeCov(X, C);

	//��������ֵ����������
	computeEig(C, vec, val);

	//������ʧ�ʣ�ȷ������ά��
	int dim = computeDim(val);
	//������
	MatrixXd res = X * vec.rightCols(dim);
	//���������ı�
	//fout << "the result is " << res.rows() << "x" << res.cols() << " after pca algorithm." << endl;
	//fout << res;


	//�����ɷַ����洢��pBuff2
	double *pBuff2 = new double[iImageSize*iBandCount];


	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < iImageSize; ++j)
		{
			pBuff2[i*iImageSize + j] = res(j, i);

		}
	}


	//�������ɷ�д��ͼ�񣨰������꼰ͶӰ��Ϣ��
	for (int i = 0; i < iBandCount; i++)
	{
		char x[]=" ";
		strcpy(x, outPath);
		char dstPath[10] = {};
		sprintf(dstPath, "%d.tif", i + 1);
		strcat(x, dstPath);
		writePcaImg(x, iWidth, iHeight, pBuff2, adfGeoTransform, prj, 7, iImageSize, i + 1);
		cout << "pca " << i + 1 << " complete" << endl;
			
	}




	//��ʾͼ��
	/*for (int i = 0; i < iBandCount; i++)
	{
	char dstPath[128] = {};
	sprintf(dstPath, "D:/pca_temp/pca%d.tif", i + 1);
	Mat mat = imread(dstPath);
	imshow("test", mat);

	}*/


	cout << "pca complete!" << endl;
	//cin.get();
	return 0;


}