# Define a range for p values from 0 to 1
p_values = np.linspace(0, 1, 100)

# Time values for t / (2N) = 1/10, 1/2, 1, 2, 10, 100
t_values = np.array([1/10, 1/2, 1, 2, 10, 100]) * (2 * N)
t_labels = ['1/10', '1/2', '1', '2', '10', '100']

# Create a plot for each t value
plt.figure(figsize=(12, 8))

for i, t in enumerate(t_values):
    # Calculate phi for each p
    phi_values = phi(p_values, t, N)
    
    # Plotting
    plt.plot(p_values, phi_values, label=f't/(2N) = {t_labels[i]}')

plt.title('Function $\phi(x, t)$ for different $t/(2N)$ values')
plt.xlabel('$p$')
plt.ylabel('$\phi(x, t)$')
plt.legend()
plt.grid(True)
plt.show()
