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
            var speed = 3.0;
            var angle = 130.2;
            var delta_t = 1;    // time units between now and then

            var filter = new DiscreteKalman2d(initial_x, initial_y, initial_cov);

            var result = filter.PredictNextLocation(speed, angle, delta_t);

            Console.WriteLine("x: " + result[0] + ", y: " + result[1]);
            Console.ReadKey();

            var measured_x = -0.5;
            var measured_y = -6.0;
            speed = 2.0;
            angle = 100;

            filter.SetLocation(measured_x, measured_y);
            result = filter.PredictNextLocation(speed, angle, delta_t);

            Console.WriteLine("x: " + result[0] + ", y: " + result[1]);
            Console.ReadKey();
        }
    }
}
