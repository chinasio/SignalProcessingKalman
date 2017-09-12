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
        /// The current measurements of the system. Y_k
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
        /// State transition matrix. 
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

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="x_location">Initial x location</param>
        /// <param name="y_location">Initial y location</param>
        public Location2d(double x_location, double y_location)
        {
            // Form initial position in column matrix
            x0 = M.DenseOfColumnMajor(2, 1, new double[] { x_location, y_location });

            // Form default initial covariance matrix in 2x2 diagonal matrix          
            SetCovarianceDiagonal(_initialCovarianceValue);

            y_k = x_k = x0;

            _logger.Info("Initial location " + x0);
        }

        /// <summary>
        /// Set the covariance matrix as a diagonal matrix
        /// </summary>
        /// <param name="value">Value on the diagonal</param>
        public void SetCovarianceDiagonal(double value)
        {
            Cov = M.Diagonal(2, 2, value);

            _logger.Info("Covariance " + Cov);
        }

        /// <summary>
        /// Set the location with the current measurement of the system
        /// </summary>
        /// <param name="x_location">location on x-axis</param>
        /// <param name="y_location">location on y-axis</param>
        public void SetMeasuredLocation(double x_location, double y_location)
        {
            y_k = Matrix<double>.Build.DenseOfColumnMajor(2, 1, new double[] { x_location, y_location });

            _logger.Info("Measured location " + Measured);
        }

        /// <summary>
        /// Set the current velocity
        /// </summary>
        /// <param name="speed">Value of the speed</param>
        /// <param name="angle">Speed direction angle</param>
        public void SetVelocity(double speed, double angle)
        {
            var x_velocity = speed * Trig.Cos(angle);
            var y_velocity = speed * Trig.Sin(angle);

            F = M.DiagonalOfDiagonalArray(new double[] { x_velocity, y_velocity });

            _logger.Info("Current velocity " + F);
        }

        /// <summary>
        /// Update measurement from the system reading.
        /// </summary>
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

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="initialLocation">Initial location object</param>
        public DiscreteKalman2d(Location2d initialLocation)
        {
            _location = initialLocation ?? throw new ArgumentException("Initial location required");
            dkf = new DiscreteKalmanFilter(initialLocation.Initial, initialLocation.Cov);
        }

        /// <summary>
        /// Constructor 
        /// </summary>
        /// <param name="initial_x">Initial x location</param>
        /// <param name="initial_y">Initial y location</param>
        /// <param name="covariance">Initial covariance value</param>
        public DiscreteKalman2d(double initial_x, double initial_y, double covariance)
        {
            _location = new Location2d(initial_x, initial_y);
            _location.SetCovarianceDiagonal(covariance);
            dkf = new DiscreteKalmanFilter(_location.Initial, _location.Cov);
        }

        /// <summary>
        /// Use the current Kalman filter and the current velocity to predict the next location  
        /// </summary>
        /// <param name="speed">System speed</param>
        /// <param name="angle">System heading angle</param>
        /// <param name="timespan">Prediction timespan</param>
        /// <returns>An array of the 2D coordinate [x, y]</returns>
        public double[] PredictNextLocation(double speed, double angle, int timespan)
        {
            // Update
            UpdateKalmanEstimate();

            // Predict
            PredictKalman(speed, angle, timespan);

            return dkf.State.ToColumnWiseArray();
        }

        /// <summary>
        /// Set the location to the current measured 
        /// </summary>
        /// <param name="x_location">x location</param>
        /// <param name="y_location">y location</param>
        public void SetLocation(double x_location, double y_location)
        {
            _location.SetMeasuredLocation(x_location, y_location);
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

        /// <summary>
        /// This step is to predict the next location based on the velocity and timespan
        /// </summary>
        /// <param name="speed">Current speed</param>
        /// <param name="angle">Current heading angle</param>
        /// <param name="timespan">How much future (discrete time units) to predict</param>
        private void PredictKalman(double speed, double angle, int timespan)
        {
            _location.SetVelocity(timespan*speed, angle);
            dkf.Predict(_location.F);
        }
    }
}
