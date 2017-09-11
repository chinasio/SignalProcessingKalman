using System;
using MathNet.Numerics;
using MathNet.Filtering.Kalman;
using MathNet.Numerics.LinearAlgebra;
using NLog;

namespace SignalProcessingKalman
{
    /// <summary>
    /// "State" of the Kalman filter
    /// </summary>
    class Location2d
    {
        private static MatrixBuilder<double> M = Matrix<double>.Build;
        private static Logger _logger = LogManager.GetCurrentClassLogger();

        private static readonly double _initialCovarianceValue = 100.0;

        private Matrix<double> x0;

        /// <summary>
        /// The initial state. X_0
        /// </summary>
        public Matrix<double> Initial
        {
            get { return x0; }
            set { x0 = value; }
        }

        private Matrix<double> x_k;

        /// <summary>
        /// The current estimated state. X_k
        /// </summary>
        public Matrix<double> Current
        {
            get { return x_k; }
            set { x_k = value; }
        }

        private Matrix<double> y_k;

        /// <summary>
        /// The measurements of the system. Y_k
        /// </summary>
        public Matrix<double> Measured
        {
            get { return y_k; }
            set { y_k = value; }
        }

        /// <summary>
        /// The covariance of the state estimate.
        /// </summary>
        public Matrix<double> Cov { get; set; }

        /// <summary>
        /// State transition matrix
        /// </summary>
        public Matrix<double> F { get; set; }

        /// <summary>
        /// Measurement model. Y_k/X_k
        /// </summary>
        public Matrix<double> H { get; set; }

        /// <summary>
        /// Covariance of measurements.
        /// </summary>
        public Matrix<double> R { get; set; }

        // Constructor
        public Location2d(double x_location, double y_location)
        {
            // Form initial position in column matrix
            x0 = M.DenseOfColumnMajor(2, 1, new double[] { x_location, y_location });

            // Form default initial covariance matrix in 2x2 diagonal matrix          
            SetCovarianceDiagonal(_initialCovarianceValue);

            y_k = x_k = x0;

            _logger.Info("Initial location " + x0);
        }

        public void SetCovarianceDiagonal(double value)
        {
            Cov = M.Diagonal(2, 2, value);

            _logger.Info("Covariance " + Cov);
        }

        public void SetMeasuredLocation(double x_location, double y_location)
        {
            y_k = Matrix<double>.Build.DenseOfColumnMajor(2, 1, new double[] { x_location, y_location });

            _logger.Info("Measured location " + Measured);
        }

        public void SetSpeed(double velocity, double angle)
        {
            var x_velocity = velocity * Trig.Cos(angle);
            var y_velocity = velocity * Trig.Sin(angle);

            F = M.DiagonalOfDiagonalArray(new double[] { x_velocity, y_velocity });

            _logger.Info("Current speed " + F);
        }

        public void UpdateMeasurement()
        {
            var s = Measured.PointwiseDivide(Current).ToColumnWiseArray();            
            H = M.Diagonal(2, 2, s);

            _logger.Info("Current measurement " + F);
        }
    }

    /// <summary>
    /// "Context" of the Kalman filter
    /// </summary>
    class DiscreteKalman2d
    {
        private DiscreteKalmanFilter dkf;
        private Location2d _location;

        public Location2d Location
        {
            get { return _location; }
            set { _location = value; }
        }

        // Constructor
        public DiscreteKalman2d(Location2d initialLocation)
        {
            _location = initialLocation ?? throw new ArgumentException("Initial location required");
            dkf = new DiscreteKalmanFilter(initialLocation.Initial, initialLocation.Cov);
        }

        public DiscreteKalman2d(double initial_x, double initial_y, double covariance)
        {
            _location = new Location2d(initial_x, initial_y);
            _location.SetCovarianceDiagonal(covariance);
            dkf = new DiscreteKalmanFilter(_location.Initial, _location.Cov);
        }

        public double[] PredictNextLocation(double velocity, double angle, int steps)
        {
            // Update
            UpdateKalmanEstimate();

            // Predict
            PredictKalman(velocity, angle, steps);

            return dkf.State.ToColumnWiseArray();
        }

        /// <summary>
        /// This step is to update the Kalman filter postior to the location measurement
        /// </summary>
        private void UpdateKalmanEstimate()
        {
            _location.R = _location.Cov.Multiply(0.3); // 0.3 is an abitrary magic number
            _location.UpdateMeasurement();

            dkf.Update(_location.Measured, _location.H, _location.R);

            _location.Current = dkf.State;
            _location.Cov = dkf.Cov;
        }

        private void PredictKalman(double velocity, double angle, int steps)
        {
            _location.SetSpeed(steps*velocity, angle);
            dkf.Predict(_location.F);
        }
    }
}
