#include "Photo.h"

void similarity_A(Matrix &A, PointList &coordinates) {
	int n = coordinates.size();
	A.resize(2 * n, 4);
	A.clear();

	for (int i = 0; i < n; i++) {
		A[2 * i][0] = coordinates[i].x;
		A[2 * i][1] = -coordinates[i].y;
		A[2 * i][2] = 1.0;

		A[2 * i + 1][0] = coordinates[i].y;
		A[2 * i + 1][1] = coordinates[i].x;
		A[2 * i + 1][3] = 1.0;
	}
}

void similarity_N(Matrix &N, PointList &coordinates) {
	int n = coordinates.size();
	N.resize(4, 4);
	N.clear();

	for (int i = 0; i < n; i++) {
		N[0][0] += coordinates[i].x * coordinates[i].x + coordinates[i].y * coordinates[i].y;
		N[0][2] += coordinates[i].x;
		N[0][3] += coordinates[i].y;
	}
	N[1][1] = N[0][0];
	N[2][2] = n;
	N[3][3] = n;

	N[2][0] = N[0][2];
	N[1][3] = N[0][2];
	N[3][1] = N[0][2];

	N[3][0] = N[0][3];
	N[1][2] = -N[0][3];
	N[2][1] = -N[0][3];
}

void similarity_u(Matrix &u, PointList &coordinates, PointList &observations) {
	int n = coordinates.size();
	u.resize(4, 1);
	u.clear();

	int xy = 0;
	for (int i = 0; i < n; i++) {
		u[0][0] += coordinates[i].x * observations[i].x + coordinates[i].y * observations[i].y;
		u[1][0] += coordinates[i].x * observations[i].y - coordinates[i].y * observations[i].x;
		u[2][0] += observations[i].x;
		u[3][0] += observations[i].y;
	}
}

void affine_A(Matrix &A, PointList &coordinates) {
	int n = coordinates.size();
	A.resize(2 * n, 6);
	A.clear();

	for (int i = 0; i < n; i++) {
		A[2 * i][0] = coordinates[i].x;
		A[2 * i][1] = coordinates[i].y;
		A[2 * i][2] = 1.0;
		A[2 * i + 1][3] = coordinates[i].x;
		A[2 * i + 1][4] = coordinates[i].y;
		A[2 * i + 1][5] = 1.0;
	}
}

void affine_N(Matrix &N, PointList &coordinates) {
	int n = coordinates.size();
	N.resize(6, 6);
	N.clear();
	
	for (int j = 0; j <= 3; j += 3) {
		for (int i = 0; i < n; i++) {
			N[j][j] += coordinates[i].x * coordinates[i].x;
			N[j][j + 1] += coordinates[i].x * coordinates[i].y;
			N[j][j + 2] += coordinates[i].x;
			N[j + 1][j + 1] += coordinates[i].y * coordinates[i].y;
			N[j + 1][j + 2] += coordinates[i].y;
		}
		N[j + 2][j + 2] = n;

		N[j + 1][j] = N[j][j + 1];
		N[j + 2][j] = N[j][j + 2];
		N[j + 2][j + 1] = N[j + 1][j + 2];
	}
}

void affine_u(Matrix &u, PointList &coordinates, PointList &observations) {
	int n = coordinates.size();
	u.resize(6, 1);
	u.clear();

	for (int i = 0; i < n; i++) {
		u[0][0] += coordinates[i].x * observations[i].x;
		u[1][0] += coordinates[i].y * observations[i].x;
		u[2][0] += observations[i].x;
		u[3][0] += coordinates[i].x * observations[i].y;
		u[4][0] += coordinates[i].y * observations[i].y;
		u[5][0] += observations[i].y;
	}
}