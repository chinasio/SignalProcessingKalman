using System;

namespace SignalProcessingKalman
{
    class Program
    {
        static void Main(string[] args)
        {
            var initial_x = 1.0;
            var initial_y = 2.0;
            var initial_cov = 4.0;
            var velocity = 3.0;
            var angel = 130.2;
            var delta_t = 1;    // time units between now and then

            var filter = new DiscreteKalman2d(initial_x, initial_y, initial_cov);

            var result = filter.PredictNextLocation(velocity, angel, delta_t);

            Console.WriteLine(result[0] + " " + result[1]);
            Console.ReadKey();
        }
    }
}
