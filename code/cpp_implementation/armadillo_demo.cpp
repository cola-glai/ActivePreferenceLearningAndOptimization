/**
 * created by gylai
 * cola lab @ uestc
 * Jan 21, 2022
 * e-mail: guiyulai AT outloock DOT com
 * 
**/

#include "base.h"

using namespace std;
using namespace arma;

void active_ranking(mat X, mat w_star);
void build_hyperplanes(mat X, cube &H);
void get_undefined_preference(Mat<int> Qh, uvec hi, int j, uvec &goodInds);
uvec get_toSort(uvec goodInds);
uvec compare_sort(uvec goodInds, Mat<int> Qhyp);
void test_compare_sort();
void quicksort_handle(Mat<int> cmp, uvec &index, int l, int r);
pair<uvec, double> linear_program_solver(mat X, Mat<int> y);
void test_linear_program_solver();
bool preference_function(vec obj1, vec obj2, vec w_star);



int main(void)
{
    uint32_t n = 10;            // number of objects
    uint32_t d = 2;             // dimensionality

    arma_rng::set_seed(0);
    // arma_rng::set_seed_random(); // set the seed to a random value

    // initialize the matrix of object
    mat X(n, d, fill::randu);

    // generate the random reference point
    mat w_star(1, d, fill::randu);

    // active_ranking(X, w_star);

    test_compare_sort();
    test_linear_program_solver();

    return 0;
}

void active_ranking(mat X, mat w_star)
{
    int n = X.n_rows; // number of objects
    int d = X.n_cols; // dimension

    printf("<info> matrix of objects(%dx%d):\n", n, d);
    cout<<X<<endl;

    printf("<info> reference point(ground truth):\n");
    cout<<w_star<<endl;


    // cube(n_rows, n_cols, n_slices, fill_form)
    // each slice is a matrix
    // n (d+1, n) matrices
    cube H(d + 1, n, n, fill::zeros);
    build_hyperplanes(X, H);


    // use random indices to randomly enumerate the objects
    // uvec hi = randperm(n);
    uvec hi = {4,1,8,7,5,2,3,0,9,6};
    cout<<"hi: "<<hi<<endl;

    // Qh: preference label of pairwise comparisons
    Mat<int> Qh = zeros<Mat<int>>(n,n);
    cout<<"Qh: "<<Qh<<endl;
    umat known;

    int jj = 0;
    while (jj < n - 1)
    {
        jj = jj + 1;

        // logical array
        uvec goodInds(jj, fill::zeros);

        // undefined preference relations related to object{hi(jj)} in Qh
        get_undefined_preference(Qh, hi, jj, goodInds);
        cout<<"goodInds: "<<goodInds<<endl;

        while(sum(goodInds) > 0)
        {
            uvec toSort = get_toSort(goodInds);
            // uvec toSort2 = linspace<uvec>(0, jj - 1, jj);
            cout<<"toSort:"<<toSort<<endl;
            Mat<int> Qhyp = Qh(hi(toSort), hi(toSort));
            cout<<"Qhyp: "<<Qhyp<<endl;


            uvec index = compare_sort(goodInds, Qhyp);
            cout<<"index: "<<index<<endl;

            uvec list = toSort(index);
            uword bis = (list.n_elem + 1) / 2;
            uword ii = list(bis);
            uvec below = list.rows(0, bis - 1);
            uvec above = list.rows(bis + 1, list.n_elem - 1);

            uword i = hi(ii);
            uword j = hi(jj);

            mat X = zeros(known.n_rows, d + 1);
            Mat<int> Y = zeros<Mat<int>>(known.n_rows, 1);
            for (uword kkk = 0; kkk < known.n_rows; ++kkk)
            {
                X.row(kkk) = (H.slice(known(kkk, 1)).col(known(kkk, 0))).t();
                Y.row(kkk) = Qh(known(kkk, 0), known(kkk, 1));
            }
            Mat<int> assumption_prefered = {1};
            pair<uvec, double> res_alpha = linear_program_solver(join_vert(X, (H.slice(j).col(i)).t()), join_vert(Y, assumption_prefered));
            Mat<int> assumption_not_prefered = {-1};
            pair<uvec, double> res_beta = linear_program_solver(join_vert(X, (H.slice(j).col(i)).t()), join_vert(Y, assumption_not_prefered));

            if (res_alpha.second < 100 && res_beta.second < 100)
            {
                cout<<"<WARNING> broke - maxiter is too low for convergence. restarting."<<endl;
                known.shed_row(known.n_elem - 1);
            }
            else if (res_alpha.second < 100)
            {// not ambiguous
                Qh(i, j) = -1;
                Qh(j, i) = -Qh(i, j);
            }
            else if (res_beta.second < 100)
            {// not ambiguous
                Qh(i, j) = 1;
                Qh(j, i) = -Qh(i, j);
            }
            else
            {// If you are here, the query is Ambiguous.
                umat new_know = {i, j};
                known = join_vert(known, new_know);
                if (preference_function(X.row(i).t(), X.row(j).t(), w_star))
                {
                    Qh(i, j) = 1;
                }
                else
                {
                    Qh(i, j) = -1;
                }

                Qh(j,i) = -Qh(i,j);

                Mat<int> assumption = {Qh(i, j)};
                pair<uvec, double> res = linear_program_solver(join_vert(X, (H.slice(j).col(i)).t()), join_vert(Y, assumption));
                uvec w = res.first;
                uvec w_approximate = w.rows(0, w.n_elem - 2) / w.row(w.n_elem - 1);
                cout<<"<info> w_approximate:"<<endl;
                cout<<w_approximate<<endl;
            }

            if (Qh(i, j) == 1)
            {
                for (uword k = 0; k < below.n_elem; ++k)
                {
                    Qh(hi(below(k)), hi(jj)) = 1;
                    Qh(hi(jj), hi(below(k))) = -1;
                }
            }
            else if(Qh(i, j) == -1)
            {
                for (uword k = 0; k < above.n_elem; ++k)
                {
                    Qh(hi(jj), hi(above(k))) = 1;
                    Qh(hi(above(k)), hi(jj )) = -1;
                }
            }

            get_undefined_preference(Qh, hi, jj, goodInds);
        }
    }
}

void build_hyperplanes(mat X, cube &H)
{
    int n = X.n_rows; // number of objects
    int d = X.n_cols; // dimension

    int records = 0;

    // build hyperplanes by pairwise points
    for (int j = 1; j < n; ++j)
    {
        for (int i = 0; i < j; ++i)
        {
            ++records;
            rowvec obj_diff(X.row(j) - X.row(i)); // row vector with size 1xd
            rowvec obj_sum(X.row(j) + X.row(i)); // row vector with size 1xd

            vec num = -0.5 * obj_sum * obj_diff.t(); // column vector with size 1x1
            vec extension(obj_diff.t()); // one column dx1
            extension.insert_rows(d, num); // at row d, insert num, d+1 rows (d+1, 1) in total; (index starts from 0)
            double normalization_value = arma::norm(arma::abs(extension));

            rowvec obj_diff_normalized = obj_diff / normalization_value;
            vec num_normalized = num / normalization_value;

            vec extension_normalized(obj_diff_normalized.t());
            extension_normalized.insert_rows(d, num_normalized);
            H.slice(j).col(i) = extension_normalized;
            cout<<H.slice(j).col(i)<<endl;
            H.slice(i).col(j) = -H.slice(j).col(i);
        }
    }

    cout<<"<info> number of hyperplanes: "<<records<<endl;
}

void get_undefined_preference(Mat<int> Qh, uvec hi, int j, uvec &goodInds)
{
    for (int i = 0; i < j; ++i)
    {
        if(Qh(hi(j), hi(i)) == 0) {
            goodInds(i) = 1;
        }
    }
}

uvec get_toSort(uvec goodInds)
{
    std::vector<uword> v;
    for(uword i = 0; i < goodInds.n_elem; ++i)
    {
        if (goodInds(i) == 1)
        {
            v.push_back(i);
        }
    }
    uvec toSort(v);
    return toSort;
}

uvec compare_sort(uvec goodInds, Mat<int> Qhyp)
{
    uvec input = linspace<uvec>(0, sum(goodInds) - 1, sum(goodInds));
    cout<<"input: "<<input<<endl;

    uvec index(input);

    Qhyp = -Qhyp;

    quicksort_handle(Qhyp, index, 0, sum(goodInds) - 1);

    return index;
}

void test_compare_sort()
{
    cout<<"##### test_compare_sort ######"<<endl;
    uvec goodInds = {1,1,1,0,0,0,1,0};
    cout<<goodInds<<endl;
    Mat<int> Qhyp = {{0,1,1,-1},
                     {-1,0,-1,-1},
                     {-1,1,0,-1},
                     {1,1,1,0}};
    cout<<Qhyp<<endl;
    uvec index = compare_sort(goodInds, Qhyp);
    cout<<"index: "<<index<<endl;
}

void quicksort_handle(Mat<int> cmp, uvec &index, int l, int r)
{
    uword i = l + 1;
    uword j = r;
    uword p = l;
    if (l >= r)
    {
        return;
    }

    while (i <= j)
    {
        switch (cmp(index(i), index(l)))
        {
            case 1:
            {
                uword tmp = index(j);
                index(j) = index(i);
                index(i) = tmp;
                j = j - 1;
                break;
            }
            case -1:
                ++i;
                break;
            default:
                ++p;
                ++i;
                break;
        }
    }

    while (i <= j)
    {
        switch (cmp(index(i), index(l)))
        {
            case 1:
            {
                uword tmp = index(j);
                index(j) = index(i);
                index(i) = tmp;
                j = j - 1;
                break;
            }
            case -1:
                ++i;
                break;
            default:
            {
                ++p;
                uword tmp = index(p);
                index(p) = index(i);
                index(i) = tmp;
                ++i;
                break;
            }
        }
    }

    // swap "less thans" with "equals"
    uword u = 0;
    uvec tmp;
    tmp.zeros(j - l + 1);
    for (uword k = p + 1; k <= j; ++k)
    {
        tmp(u) = index(k);
        ++u;
    }
    for (uword k = l; k <= p; ++k)
    {
        tmp(u) = index(k);
        ++u;
    }
    u = 0;
    for(uword k = l; k <= j; ++k)
    {
        index(k) = tmp(u);
        ++u;
    }

    quicksort_handle(cmp, index, l, l + j - p);
    quicksort_handle(cmp, index, i, r);
}

pair<uvec, double> linear_program_solver(mat X, Mat<int> Y)
{
    uword l = X.n_rows;
    uword d = X.n_cols;
    mat h = X % repmat(Y, 1, d);

    mat A = join_cols(h, -h, ones(l, 1));
    A = -A;
    mat b = ones(l, 1);
    b = -b;
    mat tmp1 = {1};
    mat f = join_vert(zeros(2 * d, 1), tmp1);
    mat tmp2 = {-1};
    mat LB = join_vert(zeros(2 * d, 1), tmp2);
    mat UB(2 *d + 1, 1);
    UB.fill(datum::inf);




}

void test_linear_program_solver()
{
    mat X = {{0.5617, 0.5006, -0.6587},
             {0.5558, 0.5616, -0.6129},
             {-0.5634, -0.3678, 0.7398},
             {-0.8803, 0.2895, 0.3758},
             {-0.2467, 0.9161, -0.3160},
             {-0.9206, 0.3607, 0.1498}};
    Mat<int> Y = ones<Mat<int>>(6, 1);
    Y(2, 0) = -1;
    Y(3, 0) = -1;
    Y(4, 0) = -1;
    cout<<"Y: "<<Y<<endl;
    pair<uvec, double> res = linear_program_solver(X, Y);
    cout<<"linear solve res: "<<res.first<<endl;
    cout<<"acc: "<<res.second<<endl;
}

bool preference_function(vec obj1, vec obj2, vec w_star)
{
    if (norm(w_star - obj1) < norm(w_star - obj2))
    {
        return true;
    }
    return false;
}