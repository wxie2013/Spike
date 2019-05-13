#include "aki_model.C"

int main()
{
    //aki_model model(false);
    aki_model model(true);

    model.run_spiking_model(false);

    return 0;
}

