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

void active_ranking(mat X, vec w_star);
void build_hyperplanes(mat X, cube &H);
void get_undefined_preference(Mat<int> Qh, uvec hi, int j, uvec &goodInds);
uvec get_toSort(uvec goodInds);
uvec compare_sort(uvec goodInds, Mat<int> Qhyp);
void test_compare_sort();
void quicksort_handle(Mat<int> cmp, uvec &index, int l, int r);
pair<vec, double> linear_program_solver(mat X, Mat<int> y);
void test_linear_program_solver();
bool preference_function(vec obj1, vec obj2, vec w_star);

// element-wise equality evaluation of two objects
bool is_equal_matrix(Mat<int> Q1, Mat<int> Q2);



int main(void)
{
//    test_compare_sort();
//    test_linear_program_solver();



    uint32_t n = 10;            // number of objects
    uint32_t d = 2;             // dimensionality

    arma_rng::set_seed(3);
    // arma_rng::set_seed_random(); // set the seed to a random value

    // initialize the matrix of object
    mat X(n, d, fill::randu);

    // generate the random reference point
    vec w_star(d, fill::randu);

    active_ranking(X, w_star);

    return 0;
}

void active_ranking(mat objects, vec w_star)
{
    vec w_approximate;
    int cnt_ambiguous_query = 0;
    int cnt_unambiguous_query = 0;
    int cnt_total_pairwise_comparison = 0;

    int n = objects.n_rows; // number of objects
    int d = objects.n_cols; // dimension

    printf("<info> matrix of objects(%dx%d):\n", n, d);
    cout<<objects<<endl;

    printf("<info> reference point(ground truth):\n");
    cout<<w_star<<endl;


    // cube(n_rows, n_cols, n_slices, fill_form)
    // each slice is a matrix
    // n (d+1, n) matrices
    cube H(d + 1, n, n, fill::zeros);
    build_hyperplanes(objects, H);


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
            cout<<"list: "<<list<<endl;
            uword bis = list.n_elem / 2;
            cout<<"bis: "<<bis<<endl;
            uword ii = list(bis);
            cout<<"ii: "<<ii<<endl;
            uvec below;
            if (bis >= 1)
            {
                below = list.rows(0, bis - 1);
            }
            uvec above;
            if (bis + 1 < list.n_elem)
            {
               above = list.rows(bis + 1, list.n_elem - 1);
            }

            cout<<"<info> get below and above finished\n";
            cout<<"below.n_elem: "<<below.n_elem<<endl;
            cout<<"below.n_rows: "<<below.n_rows<<endl;
            cout<<"above.n_rows: "<<above.n_rows<<endl;

            uword i = hi(ii);
            uword j = hi(jj);
            cnt_total_pairwise_comparison++;
            cout<<"i: "<<i<<",j: "<<j<<endl;

            mat X = zeros(known.n_rows, d + 1);
            Mat<int> Y = zeros<Mat<int>>(known.n_rows, 1);
            for (uword kkk = 0; kkk < known.n_rows; ++kkk)
            {
                X.row(kkk) = (H.slice(known(kkk, 1)).col(known(kkk, 0))).t();
                Y.row(kkk) = Qh(known(kkk, 0), known(kkk, 1));
            }
            Mat<int> assumption_prefered = {1};
            pair<vec, double> res_alpha = linear_program_solver(join_vert(X, (H.slice(j).col(i)).t()), join_vert(Y, assumption_prefered));
            Mat<int> assumption_not_prefered = {-1};
            pair<vec, double> res_beta = linear_program_solver(join_vert(X, (H.slice(j).col(i)).t()), join_vert(Y, assumption_not_prefered));

            cout<<"^"<<endl;
            cout<<"res_alpha.second: "<<res_alpha.second<<endl;
            cout<<"res_beta.second: "<<res_beta.second<<endl;
            if (res_alpha.second < 100 && res_beta.second < 100)
            {
                cout<<"<WARNING> broke - maxiter is too low for convergence. restarting."<<endl;
                known.shed_row(known.n_elem - 1);
            }
            else if (res_alpha.second < 100)
            {// not ambiguous
                Qh(i, j) = -1;
                Qh(j, i) = -Qh(i, j);
                cnt_unambiguous_query++;
            }
            else if (res_beta.second < 100)
            {// not ambiguous
                Qh(i, j) = 1;
                Qh(j, i) = -Qh(i, j);
                cnt_unambiguous_query++;
            }
            else
            {// If you are here, the query is Ambiguous.
                cnt_ambiguous_query++;
                cout<<"<info> the query is Ambiguous"<<endl;
                umat new_know = {i, j};
                known = join_vert(known, new_know);
                cout<<"known: "<<known<<endl;
                cout<<objects.row(i)<<endl;
                if (preference_function(objects.row(i).t(), objects.row(j).t(), w_star))
                {
                    Qh(i, j) = 1;
                }
                else
                {
                    Qh(i, j) = -1;
                }

                Qh(j,i) = -Qh(i,j);

                Mat<int> assumption = {Qh(i, j)};
                pair<vec, double> res = linear_program_solver(join_vert(X, (H.slice(j).col(i)).t()), join_vert(Y, assumption));
                cout<<"<info> w_approximate:"<<endl;
                w_approximate = res.first;
                cout<<res.first<<endl;
            }

            if (Qh(i, j) == 1)
            {
                cout<<"<info> induced by below quick sort: "<<below.n_rows<<endl;
                cnt_unambiguous_query += below.n_rows;
                cnt_total_pairwise_comparison += below.n_rows;
                for (uword k = 0; k < below.n_rows; ++k)
                {
                    if (hi(below(k)) == hi(jj))
                    {
                        continue;
                    }
                    Qh(hi(below(k)), hi(jj)) = 1;
                    Qh(hi(jj), hi(below(k))) = -1;
                }
            }
            else if(Qh(i, j) == -1)
            {
                cout<<"<info> induced by above quick sort: "<<above.n_rows<<endl;
                cnt_unambiguous_query += above.n_rows;
                cnt_total_pairwise_comparison += above.n_rows;
                for (uword k = 0; k < above.n_rows; ++k)
                {
                    if (hi(above(k)) == hi(jj))
                    {
                        continue;
                    }
                    Qh(hi(jj), hi(above(k))) = 1;
                    Qh(hi(above(k)), hi(jj )) = -1;
                }
            }

            get_undefined_preference(Qh, hi, jj, goodInds);
            cout<<"goodInds: "<<goodInds<<endl;
        }
    }

    cout<<"cnt_ambiguous_query: "<<cnt_ambiguous_query<<endl;
    cout<<"cnt_unambiguous_query: "<<cnt_unambiguous_query<<endl;
    cout<<"cnt_total_pairwise_comparison: "<<cnt_total_pairwise_comparison<<endl;
    cout<<"w_approximate: "<<w_approximate<<endl;

    cout<<"Qh after active ranking:"<<endl;
    cout<<Qh<<endl;

    Mat<int> Q_gold = zeros<Mat<int>>(n,n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (i == j)
            {
                continue;
            }
            if (preference_function(objects.row(i).t(), objects.row(j).t(), w_approximate))
            {
                Q_gold(i, j) = 1;
                Q_gold(j, i) = -1;
            }
            else
            {
                Q_gold(i, j) = -1;
                Q_gold(j, i) = 1;
            }
        }
    }
    cout<<"Q_gold:"<<endl;
    cout<<Q_gold<<endl;
    if (is_equal_matrix(Qh, Q_gold))
    {
        cout<<"<info> preference matrix induced by w_star is equal to which induced by w_approximate!"<<endl;
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
        else
        {
            goodInds(i) = 0;
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

pair<vec, double> linear_program_solver(mat X, Mat<int> Y)
{
    uword l = X.n_rows;
    uword d = X.n_cols;
    mat h = X % repmat(Y, 1, d);
    cout<<"size of h: "<<h.n_rows<<","<<h.n_cols<<endl;

    mat b = ones(l, 1);
    cout<<"size of b: "<<b.n_rows<<","<<b.n_cols<<endl;
    mat A = join_rows(h, -h, b);
    A = -A;
    b = -b;
    mat tmp1 = {1};
    cout<<"join_cols ok\n";
    mat f = join_vert(zeros(2 * d, 1), tmp1);
    mat tmp2 = {-1};
    mat LB = join_vert(zeros(2 * d, 1), tmp2);
    mat UB(2 *d + 1, 1);
    UB.fill(datum::inf);

    // lp_solve solves the linear programming solver based on the revised simplex method
    lprec *lp;
//    int j;
    int Ncol, *colno = NULL, ret = 0;
    REAL *row = NULL;

    Ncol = 2 * d + 1; // number of variables
    lp = make_lp(0, Ncol);

    if (lp == NULL)
    {
        ret = 1; /* couldn't construct a new model... */
    }

    if (ret == 0)
    {
        /* let us name our variables. Not required, but can be useful for debugging */
        // set_col_name(lp, 1, "x");
        // set_col_name(lp, 2, "y");

        /* create space large enough for one row */
        colno = (int *)malloc(Ncol * sizeof(*colno));
        row = (REAL *)malloc(Ncol * sizeof(*row));
        if ((colno == NULL) || (row == NULL))
        {
            ret = 2;
        }
    }

    if (ret == 0)
    {
        set_add_rowmode(lp, TRUE);
        // add constraint: A*u <= b
        for (uword i = 0; i < l; ++i)
        {
            int j = 0;
            for (; j < Ncol; ++j)
            {
                colno[j] = j + 1;
                row[j] = A(i, j);
            }
            /* add the row to lpsolve */
            if (!add_constraintex(lp, j, row, colno, LE, b(i, 0)))
            {
                ret = 3;
                break;
            }
        }

//        /* construct first row (120 x + 210 y <= 15000) */
//        j = 0;
//        colno[j] = 1; /* first column */
//        row[j++] = 120;
//        colno[j] = 2; /* second column */
//        row[j++] = 210;
//
//        /* add the row to lpsolve */
//        if (!add_constraintex(lp, j, row, colno, LE, 15000))
//        {
//            ret = 3;
//        }
    }

    if (ret == 0)
    {
        // add constraint: lower bound [0,0,0,...,-1]
        mat M = eye(Ncol, Ncol);
        mat lb = zeros(Ncol, 1);
        lb(lb.n_rows - 1, 0) = -1;
        for (int i = 0; i < Ncol; ++i)
        {
            int k = 0;
            for (; k < Ncol; ++k)
            {
                colno[k] = k + 1;
                row[k] = M(i, k);
            }
            /* add the row to lpsolve */
            if (!add_constraintex(lp, k, row, colno, GE, lb(i, 0)))
            {
                ret = 3;
                break;
            }
        }
//        /* construct second row (110 x + 30 y <= 4000) */
//        j = 0;
//        colno[j] = 1; /* first column */
//        row[j++] = 110;
//        colno[j] = 2; /* second column */
//        row[j++] = 30;
//        /* add the row to lpsolve */
//        if (!add_constraintex(lp, j, row, colno, LE, 4000))
//        {
//            ret = 3;
//        }
    }

//    if (ret == 0)
//    {
//        /* construct third row (x + y <= 75) */
//        j = 0;
//        colno[j] = 1; /* first column */
//        row[j++] = 1;
//        colno[j] = 2; /* second column */
//        row[j++] = 1;
//        /* add the row to lpsolve */
//        if (!add_constraintex(lp, j, row, colno, LE, 75))
//        {
//            ret = 3;
//        }
//    }
    if (ret == 0)
    {
        set_add_rowmode(lp, FALSE);
        /* set the objective function (u_(2*d + 1)) */
        int i = 0;
        for (; i < Ncol; ++i)
        {
            colno[i] = i + 1;
            row[i] = f(i, 0);
        }
        if (!set_obj_fnex(lp, i, row, colno))
        {
            ret = 4;
        }
//        j = 0;
//        colno[j] = 1; /* first column */
//        row[j++] = 143;
//        colno[j] = 2; /* second column */
//        row[j++] = 60;
//        if (!set_obj_fnex(lp, j, row, colno))
//        {
//            ret = 4;
//        }
    }
    if (ret == 0)
    {
        /* set the object direction to minimize */
        set_minim(lp);
        write_LP(lp, stdout);
        set_verbose(lp, IMPORTANT);
        ret = solve(lp);
        if (ret == OPTIMAL)
        {
            cout<<"<info> get optimal solution in linear programming solving\n";
            ret = 0;
        }
        else
        {
            ret = 5;
        }
    }
    if (ret == 0)
    {
        /* objective value */
        printf("Objective value: %f", get_objective(lp));
        /* variable values */
        get_variables(lp, row);
        for (int j = 0; j < Ncol; j++)
        {
            printf("%s: %f", get_col_name(lp, j + 1), row[j]);
        }
        cout<<endl;
    }

    if (colno != NULL)
    {
        free(colno);
        colno = NULL;
    }
    if (lp != NULL)
    {
        delete_lp(lp);
        lp = NULL;
    }

    vec u1 = zeros(d);
    vec u2 = zeros(d);
    for (uword i = 0; i < d; ++i)
    {
        u1(i) = row[i];
    }
    for (uword i = 0; i < d; ++i)
    {
        u2(i) = row[d + i];
    }
    vec w = u1 - u2;
    vec s = X * w;
    int number_right_judges = 0;
    for (uword i = 0; i < l; ++i)
    {
        if ((Y(i, 0) == -1 && s(i) < 0) ||
            (Y(i, 0) == 1 && s(i) > 0))
        {
            ++number_right_judges;
        }
    }

    double acc = (double)number_right_judges / l * 100;
    cout<<"acc: "<<acc<<endl;

    vec w_approximate;
    if (fabs(w(w.n_rows - 1)) < 1e-6)
    {// handle <divide zero>
        cout<<"<WARNING> divide zero"<<endl;
        w_approximate = zeros(d - 1, 1);
    }
    else
    {
        w_approximate = w.rows(0, w.n_rows - 2) / w(w.n_rows - 1);
    }

    pair<vec, double> res(w_approximate, acc);

    if (row != NULL)
    {
        free(row);
        row = NULL;
    }

    return res;
}

void test_linear_program_solver()
{
    mat X = {{0.5617, 0.5006, -0.6587},
             {0.5558, 0.5616, -0.6129},
             {-0.5634, -0.3678, 0.7398},
             {-0.8803, 0.2895, 0.3758},
             {-0.2467, 0.9161, -0.3160},
             {-0.9206, 0.3607, 0.1498},
             {-0.7321, -0.5009, 0.4616},
             {-0.5714, 0.7682, -0.2889},
             {-0.0299, 0.8243, -0.5653},
             {-0.8318, 0.5091, -0.2212},
             {0.6702, -0.7039, 0.2352},
             {0.4582, 0.7453, -0.4843},
             {-0.5954, -0.6146, 0.5174},
             {0.9784, -0.2060, -0.0168},
             {-0.1882, -0.8074, 0.5593},
             {0.0649, 0.8503, -0.5222},
             {0.1401, 0.8415, -0.5218}
             };

    Mat<int> Y = ones<Mat<int>>(17, 1);
    Y(2, 0) = -1;
    Y(3, 0) = -1;
    Y(4, 0) = -1;
    Y(5, 0) = -1;
    Y(6, 0) = -1;
    Y(7, 0) = -1;
    Y(11, 0) = -1;
    Y(12, 0) = -1;
    Y(13, 0) = -1;
    Y(14, 0) = -1;

    linear_program_solver(X, Y);
}

bool preference_function(vec obj1, vec obj2, vec w_star)
{
    if (norm(w_star - obj1) < norm(w_star - obj2))
    {
        return true;
    }
    return false;
}

bool is_equal_matrix(Mat<int> Q1, Mat<int> Q2)
{
    umat is_equal = Q1 == Q2;
    cout<<"is_equal:"<<endl;
    cout<<is_equal<<endl;
    uvec non_zero_row = all(is_equal, 1);
    return all(non_zero_row);
}