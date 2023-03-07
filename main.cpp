// calculating distance field
tira::image<float> dist(tira::image<float> levelset) {


	tira::image<float> binary_boundary(levelset.width(), levelset.height());

	std::vector<std::tuple<int, int>> neighbors;
	neighbors.emplace_back(0, 1);
	neighbors.emplace_back(1, 0);
	neighbors.emplace_back(-1, 0);
	neighbors.emplace_back(0, -1);


	// indentifying boundary cells
	for (int y = 1; y < binary_boundary.height() - 1; y++){						// for every row in the image
		for (int x = 1; x < binary_boundary.width() - 1; x++){					// for every column in the image
			for (int k=0; k < neighbors.size(); k++){				// for every neighbor

				int nx = x + get<0>(neighbors[k]);					// calculate the x coordinate of the neighbor cell
				int ny = y + get<1>(neighbors[k]);					// calculate the y coordinate of the neighbor cell

				if (levelset(x, y) * levelset(nx, ny) <= 0) {					// if the product of the current cell and neighboring cell is negative
					binary_boundary(x, y) = 1;								// this cell is a boundary cell
					binary_boundary(nx, ny) = 1;							// the neighboring cell is a boundary cell
				}				
			}
		}
	}



	tira::image<float> binary_boundary_new(levelset.width(), levelset.height());

	//// multiplying df_dist with large values
	for (int y = 0; y < binary_boundary.height(); y++)
	{
		for (int x = 0; x < binary_boundary.width(); x++)
		{

			binary_boundary_new(x, y) = 9999 + binary_boundary(x, y);
			
		}
	}



	tira::image<float> boundary(levelset.width(), levelset.height());



	// calculate the distance from boundary cells to the contour
	for (int y = 1; y < binary_boundary.height() - 1; y++) {						// for every row in the image
		for (int x = 1; x < binary_boundary.width() - 1; x++) {						// for every column in the image
			if (binary_boundary(x,y) == 1) {									// if the pixel (x, y) is in the boundary
			    for (int k = 0; k < neighbors.size(); k++) {				// for every neighbor

				   int nx = x + get<0>(neighbors[k]);						// calculate the x coordinate of the neighbor cell
				   int ny = y + get<1>(neighbors[k]);						// calculate the y coordinate of the neighbor cell
				   if (binary_boundary(nx, ny) == 1) {
					   //if (phi(x, y) * phi(nx, ny) <= 0) {
					   //if (df_dist(x, y) * df_dist(nx, ny) <= 0) {					// if the product of the current cell and neighboring cell is negative
						   float da = (abs(levelset(x, y))) / (abs(levelset(nx, ny) - levelset(x, y)));
						   float db =( abs(levelset(nx, ny))) / (abs(levelset(nx, ny) - levelset(x, y)));
						   binary_boundary_new(x, y) = std::min(binary_boundary_new(x, y), da);
						   binary_boundary_new(nx, ny) = std::min(binary_boundary_new(nx, ny), db);
					  // }
					   
				   }
				}
			}
		
		}
	}

	std::vector<float>distGrid;
	distGrid.resize(levelset.height() * levelset.width());

	const int height = levelset.height();
	const int width = levelset.width();
	const int row = width;

	for (int y = 0; y < binary_boundary_new.height(); y++)
	{
		for (int x = 0; x < binary_boundary_new.width(); x++)
		{
			distGrid[y * width + x] = binary_boundary_new(x, y);
		}

	}

	std::vector<float>frozenCells;
	frozenCells.resize(levelset.height() * levelset.width());

	for (int i = 0; i < distGrid.size(); i++)
	{
		if (distGrid[i] < 1) {
			frozenCells[i] = true;
		}
		else {
			frozenCells[i] = false;
		}

	}


	const int NSweeps = 4;



	//// sweep directions { start, end, step }
	const int dirX[NSweeps][3] = { {0, width - 1, 1} , {width - 1, 0, -1}, {width - 1, 0, -1}, {0, width - 1, 1} };
	const int dirY[NSweeps][3] = { {0, height - 1, 1}, {0, height - 1, 1}, {height - 1, 0, -1}, {height - 1, 0, -1} };
	double aa[2];
	double d_new, a, b;
	int s, ix, iy, gridPos;
	const double h = 1.0, f = 1.0;

	for (s = 0; s < NSweeps; s++) {

		for (iy = dirY[s][0]; dirY[s][2] * iy <= dirY[s][1]; iy += dirY[s][2]) {
			for (ix = dirX[s][0]; dirX[s][2] * ix <= dirX[s][1]; ix += dirX[s][2]) {

				gridPos = iy * row + ix;

				if (!frozenCells[gridPos]) {

					// === neighboring cells (Upwind Godunov) ===


					if (iy == 0 || iy == (height - 1)) {                    // calculation for ymin
						if (iy == 0) {
							aa[1] = distGrid[(iy + 1) * row + ix];
						}
						if (iy == (height - 1)) {
							aa[1] = distGrid[(iy - 1) * row + ix];
						}
					}
					else {
						aa[1] = distGrid[(iy - 1) * row + ix] < distGrid[(iy + 1) * row + ix] ? distGrid[(iy - 1) * row + ix] : distGrid[(iy + 1) * row + ix];
						//aa[1] = std::min(distGrid[(iy - 1) * row + ix], distGrid[(iy + 1) * row + ix]);
					}

					if (ix == 0 || ix == (width - 1)) {                    // calculation for xmin
						if (ix == 0) {
							aa[0] = distGrid[iy * row + (ix + 1)];
						}
						if (ix == (width - 1)) {
							aa[0] = distGrid[iy * row + (ix - 1)];
						}
					}
					else {
						aa[0] = distGrid[iy * row + (ix - 1)] < distGrid[iy * row + (ix + 1)] ? distGrid[iy * row + (ix - 1)] : distGrid[iy * row + (ix + 1)];
						//aa[0] = std::min(distGrid[iy * row + (ix - 1)], distGrid[iy * row + (ix + 1)]);
					}

					a = aa[0]; b = aa[1];
					d_new = (fabs(a - b) < f * h ? (a + b + sqrt(2.0 * f * f * h * h - (a - b) * (a - b))) * 0.5 : std::fminf(a, b) + f * h);

					distGrid[gridPos] = d_new;
				}
			}
		}
	}

	for (int y = 0; y < levelset.height(); y++)
	{
		for (int x = 0; x < levelset.width(); x++)
		{
			levelset(x, y) = distGrid[y * width + x];
		}

	}

	return levelset;
}




//calculating signed distance field
tira::image<float> calc_sdf(tira::image<float> distance) {
	//tira::image<float> SDF = distance;

	std::vector<float>SDF;
	SDF.resize(distance.height() * distance.width());

	int width = distance.width();
	int height = distance.height();

	for (int y = 0; y < distance.height(); y++)
	{
		for (int x = 0; x < distance.width(); x++)
		{
			SDF[y * width + x] = distance(x, y);
		}

	}

	std::vector <float> frozenCells;
	frozenCells.resize(distance.height() * distance.width());

	for (int i = 0; i < frozenCells.size(); i++)
	{
		if (SDF[i] < 0.5) {
			frozenCells[i] = true;
		}
		else {
			frozenCells[i] = false;
		}

	}

	for (int i = 0; i < SDF.size(); i++) {
		SDF[i] = -1 * SDF[i];
	}

	//cout << "STARTED \n";
	double val; int gridPos;
	const int row = distance.width();
	const int nx = distance.width() - 1;
	const int ny = distance.height() - 1;
	int ix = 0, iy = 0;

	std::stack<std::tuple<int, int>> stack = {};

	std::tuple<int, int> idsPair;
	// find the first unfrozen cell
	gridPos = 0;
	/*bool found = false;
	for (int ix = 0; ix < nx; ix++) {
		for (int iy = 0; iy < ny; iy++) {
			if (distGrid(ix, iy) == 0) {
				cout << "FOUND \n";
				found = true;
				break;
			}
		}
		if (found == true) {
			break;
		}
	}*/

	while (frozenCells[gridPos]) {
		ix += (ix < nx ? 1 : 0);
		iy += (iy < ny ? 1 : 0);
		gridPos = row * iy + ix;
	}
	stack.push({ ix, iy });
	// a simple pixel flood
	while (stack.size()) {
		idsPair = stack.top();
		stack.pop();
		ix = std::get<0>(idsPair);
		iy = std::get<1>(idsPair);
		gridPos = row * iy + ix;
		if (!frozenCells[gridPos]) {
			val = -1.0 * SDF[gridPos];
			SDF[gridPos] = val;
			frozenCells[gridPos] = true; // freeze cell when done
			if (ix > 0) {
				stack.push({ ix - 1, iy });
			}
			if (ix < nx) {
				stack.push({ ix + 1, iy });
			}
			if (iy > 0) {
				stack.push({ ix, iy - 1 });
			}
			if (iy < ny) {
				stack.push({ ix, iy + 1 });
			}
		}
	}


	for (int y = 0; y < distance.height(); y++)
	{
		for (int x = 0; x < distance.width(); x++)
		{
			distance(x, y) = SDF[y * width + x];
		}

	}
	return distance;
}
