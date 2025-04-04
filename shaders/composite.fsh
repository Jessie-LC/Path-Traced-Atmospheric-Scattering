#version 450

/*
    Copyright (C) 2023-2025  Jessie

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#define CIE_VERSION 0 //[0 1]
#define MAXIMUM_SAMPLE_COUNT -1 //[-1 100 200 300 400 500 1000 2000 3000 4000 5000 7000 8000 9000 10000 15000 20000 25000 30000 35000 40000 45000 50000 55000 60000 65000 70000 75000 80000 85000 90000 95000 100000] Set to -1 for unlimited sample count

//#define USE_RAYMARCHED_ATMOSPHERE

//#define ENABLE_CLOUDS //Very slow

#define RAYLEIGH_DENSITY_PROFILE 0 //[0 1]
#define AEROSOL_DENSITY_PROFILE 0 //[0 1 2]
#define OZONE_DENSITY_PROFILE 0 //[0 1 2]

#define PHASE_FUNCTION_RAYLEIGH 0 //[0 1]
#define PHASE_FUNCTION_AEROSOL 0 //[0 1 2 3 4]
#define PHASE_FUNCTION_CLOUD 0 //[0 1]

#define TURBIDITY 1 //[1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6 6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8 6.9 7 7.1 7.2 7.3 7.4 7.5 7.6 7.7 7.7 7.8 7.9 8 8.2 8.3 8.4 8.5 8.6 8.7 8.8 8.9 9 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9 10 10.1 10.2 10.3 10.4 10.5 10.6 10.7 10.8 10.9 11 11.1 11.2 11.3 11.4 11.5 11.6 11.7 11.8 11.9 12 12.1 12.2 12.3 12.4 12.5 12.6 12.7 12.8 12.9 13 13.1 13.2 13.3 13.4 13.5 13.6 13.7 13.8 13.9 14 14.1 14.2 14.3 14.4 14.5 14.6 14.7 14.8 14.9 15]

#define SUN_ANGLE 4 //[-6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120]

#define SUN_ROTATION 0 //[-10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120]

//#define VIEW_FROM_SPACE

#define PROJECTION 2 //[0 1 2]

#define SCATTERING_EVENTS 10 //[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200]

#define GROUND_ALBEDO_R 79 //[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255]
#define GROUND_ALBEDO_G 79 //[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255]
#define GROUND_ALBEDO_B 79 //[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255]

#define OZONE_MODE 0 //[0 1 2]
#define OZONE_DENSITY_MULTIPLIER 1.0 //[0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8 6.9 7.0 7.1 7.2 7.3 7.4 7.5 7.6 7.7 7.8 7.9 8.0 8.1 8.2 8.3 8.4 8.5 8.6 8.7 8.8 8.9 9.0 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9 10.0]

//#define USE_MEASURED_AEROSOL_COEFFICIENT
//#define USER_DEFINED_COEFFICIENTS

#define RAYLEIGH_COLOR_R 44 //[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255]
#define RAYLEIGH_COLOR_G 105 //[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255]
#define RAYLEIGH_COLOR_B 255 //[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255]

#define AEROSOL_COLOR_R 171 //[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255]
#define AEROSOL_COLOR_G 206 //[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255]
#define AEROSOL_COLOR_B 255 //[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255]

#define OZONE_COLOR_R 255 //[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255]
#define OZONE_COLOR_G 147 //[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255]
#define OZONE_COLOR_B 10 //[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255]

#define FOV 60 //[10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90]

const float sunRadius        = 6.95e8;
const float sunDistance      = 1.49e11;
const float sunAngularRadius = sunRadius / sunDistance;

layout(location = 0) out vec4 simulationOutput;

layout(location = 0) in vec2 textureCoordinate;

uniform sampler2D colortex0;
uniform sampler2D noisetex;

uniform sampler2D cdfTextureAerosol;
uniform sampler2D cdfTextureAerosolLowAltitude;
uniform sampler2D cdfTextureRainbow;
uniform sampler2D cdfTextureCloud;
uniform sampler2D cdfTextureRayleigh;
uniform sampler2D phaseTextureAerosol;
uniform sampler2D phaseTextureAerosolLowAltitude;
uniform sampler2D phaseTextureRainbow;
uniform sampler2D phaseTextureCloud;
uniform sampler2D phaseTextureRayleigh;
uniform sampler2D earthDiffuse;

uniform sampler2D usStandardAtmosphere;

uniform sampler3D noise3D;

uniform mat4 gbufferModelViewInverse;
uniform mat4 gbufferProjection;
uniform mat4 gbufferPreviousProjection;
uniform mat4 gbufferModelView;
uniform mat4 gbufferPreviousModelView;

uniform vec3 cameraPosition, previousCameraPosition;

uniform float viewHeight, viewWidth;

uniform int frameCounter;

const bool colortex0Clear  = false;

#include "/lib/universal/universal.glsl"
#include "/lib/rng/pcg.glsl"

#include "/lib/space/sun.glsl"
#include "/lib/atmosphere/phase.glsl"
#include "/lib/atmosphere/constants.glsl"
#include "/lib/atmosphere/misc.glsl"
#include "/lib/atmosphere/pt.glsl"
//#include "/lib/atmosphere/pt_plane.glsl"
#include "/lib/atmosphere/raymarched.glsl"

const float sensorWidth = 1e-3 * 36.0;
const float sensorHeight = 1e-3 * 27.0;

float BinarySearch_CIE2012(float lowBound, float highBound, float toFind, in vec3 weights) {
    float midIndex = (lowBound + highBound) / 2.0;

    int lookupRes = textureSize(CIELUT, 0).x * 2;
    for (int x = 0; x < log2(lookupRes); ++x) {
        float coordinate = midIndex;
        float coordTweaked = (coordinate * float(lookupRes - 1) + 0.5) / lookupRes;

        float value = saturate(weights.x * texture(CIELUT, vec2(coordTweaked, 0)).x + weights.y * texture(CIELUT, vec2(coordTweaked, 0)).y + weights.z * texture(CIELUT, vec2(coordTweaked, 0)).z); //The CDF value should never go above 1.

        if (value < toFind) {
            lowBound = midIndex;
        }
        else if (value > toFind) {
            highBound = midIndex;
        }
        else {
            return midIndex;
        }

        midIndex = (lowBound + highBound) / 2.0;
    }
    return midIndex;
}
float BinarySearch_CIE1931(float lowBound, float highBound, float toFind, in vec3 weights) {
    float midIndex = (lowBound + highBound) / 2.0;

    int lookupRes = textureSize(CIELUT_1931, 0).x * 2;
    for (int x = 0; x < log2(lookupRes); ++x) {
        float coordinate = midIndex;
        float coordTweaked = (coordinate * float(lookupRes - 1) + 0.5) / lookupRes;

        float value = saturate(weights.x * texture(CIELUT_1931, vec2(coordTweaked, 0)).x + weights.y * texture(CIELUT_1931, vec2(coordTweaked, 0)).y + weights.z * texture(CIELUT_1931, vec2(coordTweaked, 0)).z); //The CDF value should never go above 1.

        if (value < toFind) {
            lowBound = midIndex;
        }
        else if (value > toFind) {
            highBound = midIndex;
        }
        else {
            return midIndex;
        }

        midIndex = (lowBound + highBound) / 2.0;
    }
    return midIndex;
}

vec3 FisheyeProjection(in float cx, in float cy, in float fov) {
    vec2 coordinate = vec2(cx, cy);
    vec2 aspect_ratio = vec2(viewWidth / viewHeight, 1);
	vec3 p = vec3((2.0 * coordinate - 1.0) * aspect_ratio * fov, -1.0);

	float z2 = p.x * p.x + p.y * p.y;

    if (z2>1.0) {
        return vec3(-1.0);
    }

	float phi = atan(p.y, p.x);
	float theta = acos(1.0 - z2);

    vec3 dir = vec3(sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi));

    return dir;
}

vec3 EquirectangularProjection(in vec2 coord) {
    coord.x -= 0.25;
	const vec2 coordToLongLat = vec2(2.0 * pi, pi);
	coord.y -= 0.5;
	vec2 longLat = coord * coordToLongLat;
	float longitude = longLat.x;
	float latitude = longLat.y - (2.0 * pi);

	float cosLat = cos(latitude);
	float cosLong = cos(longitude);
	float sinLat = sin(latitude);
	float sinLong = sin(longitude);

	return normalize(vec3(cosLat * cosLong, cosLat * sinLong, sinLat).xzy);
}

void main() {
    int samples = int(texture(colortex0, textureCoordinate).a);

    #ifdef VIEW_FROM_SPACE
        bool    moved  = gbufferProjection != gbufferPreviousProjection;
        moved = moved || gbufferModelView  != gbufferPreviousModelView;
        moved = moved || cameraPosition    != previousCameraPosition;
        
        if(moved) {
            samples = 0;
        }
    #endif
    #ifndef VIEW_FROM_SPACE
        #if PROJECTION == 0
            bool    moved  = gbufferProjection != gbufferPreviousProjection;
            moved = moved || gbufferModelView  != gbufferPreviousModelView;
            moved = moved || cameraPosition    != previousCameraPosition;
            
            if(moved) {
                samples = 0;
            }
        #endif
    #endif
    
    #if MAXIMUM_SAMPLE_COUNT != -1
        if(samples >= MAXIMUM_SAMPLE_COUNT) {
            samples = MAXIMUM_SAMPLE_COUNT;
            vec3 previousColor = texture(colortex0, textureCoordinate).rgb;
            simulationOutput.rgb = previousColor;
            simulationOutput.a   = samples;
            return;
        }
    #endif

    const int TILE_SIZE = 540;
    const int TILE_SAMPLES = 5;

    //int horWrap = int(ceil(viewWidth / float(TILE_SIZE)));
    //int vertWrap = int(ceil(viewHeight / float(TILE_SIZE)));
    //int tile = (frameCounter / TILE_SAMPLES) % (horWrap * vertWrap);
    //ivec2 tilePosition = TILE_SIZE * ivec2(tile % horWrap, tile / horWrap);
    //bool inTile = ((gl_FragCoord.x - tilePosition.x) - TILE_SIZE / 2 > -TILE_SIZE / 2  &&
    //               (gl_FragCoord.x - tilePosition.x) - TILE_SIZE / 2 <  TILE_SIZE / 2) &&
    //              ((gl_FragCoord.y - tilePosition.y) - TILE_SIZE / 2 > -TILE_SIZE / 2  &&
    //               (gl_FragCoord.y - tilePosition.y) - TILE_SIZE / 2 <  TILE_SIZE / 2);
    //if (!inTile) {
    //    simulationOutput.rgb = texture(colortex0, textureCoordinate).rgb;
    //    simulationOutput.a   = texture(colortex0, textureCoordinate).a;
    //    return;
    //}

	uint seed = uint(gl_FragCoord.x * viewHeight + gl_FragCoord.y);
	     seed = seed * 720720u + uint(samples);

    InitRand(seed);

    vec3 sunVector = Rotate(vec3(0.0, sin(radians(SUN_ANGLE)), cos(radians(SUN_ANGLE))), vec3(0.0, 1.0, 0.0), radians(SUN_ROTATION));

    vec2 viewResolution = vec2(viewWidth, viewHeight);

    vec2 aa = (2.0 * RandNext2F() - 1.0) * rcp(viewResolution.xy);

    vec2 sensorSize = vec2(sensorWidth, sensorHeight);

    float aspectRatio = viewResolution.x / viewResolution.y;

    vec2 focalLength  = 0.5 * sensorSize / (2.0 * vec2(aspectRatio, 1.0) * tan(radians(FOV) / 2.0));

    vec3 history = vec3(0.0);
    int weight = 0;
    for(int x = -2; x <= 2; ++x) {
        for(int y = -2; y <= 2; ++y) {
            history += texture(colortex0, textureCoordinate + vec2(x, y) / viewResolution, 0).rgb;
        }
    }
    history = max(history, 0.0);
    vec3 cmfWeights = mix(saturate(history / dot(vec3(1.0), history)), vec3(1.0 / 3.0), history == vec3(0.0) ? 1.0 : 0.2);
    #if CIE_VERSION == 0
        float wavelength = (441.0 * BinarySearch_CIE2012(0.0, 1.0, RandNextF(), cmfWeights)) + 390.0;
    #elif CIE_VERSION == 1
        float wavelength = (471.0 * BinarySearch_CIE1931(0.0, 1.0, RandNextF(), cmfWeights)) + 360.0;
    #endif

    #if CIE_VERSION == 0
        vec3 cie = texture(CIELUT, vec2(int(wavelength - 390.0) / 441.0, 1)).xyz;
        float wavelengthPDF = cmfWeights.x * (cie.x / 113.042) + cmfWeights.y * (cie.y / 113.042) + cmfWeights.z * (cie.z / 113.042);
    #elif CIE_VERSION == 1
        vec3 cie = texture(CIELUT_1931, vec2(int(wavelength - 360.0) / 471.0, 1)).xyz;
        float wavelengthPDF = cmfWeights.x * (cie.x / 106.857) + cmfWeights.y * (cie.y / 106.857) + cmfWeights.z * (cie.z / 106.857);
    #endif

    #ifndef VIEW_FROM_SPACE
        vec3 viewPosition = vec3(0.0, planetRadius + 1.0, 0.0);

        #if PROJECTION == 0
            vec2 uv  = (gl_FragCoord.xy * rcp(viewResolution.xy)) * 2.0 - 1.0;
                 uv += aa;;
                 uv *= sensorSize / (2.0 * focalLength);

            vec3 viewDirection = normalize(mat3(gbufferModelViewInverse) * vec3(uv.x, uv.y, -1.0));
        #elif PROJECTION == 1
            vec2 uv  = gl_FragCoord.xy * rcp(viewResolution.xy);
                 uv += aa;

            vec3 viewDirection = EquirectangularProjection(uv);
        #else
            vec2 uv  = gl_FragCoord.xy * rcp(viewResolution.xy);
                 uv += aa;

            float fov = tan(radians(45.0));
            vec3 viewDirection = FisheyeProjection(uv.x, uv.y, fov);

            if(viewDirection == vec3(-1.0)) {
                discard;
            }
        #endif
    #else
        vec3 viewPosition = vec3(0.0, 0.0, planetRadius * 4.0);

        vec2 uv  = (gl_FragCoord.xy * rcp(viewResolution.xy)) * 2.0 - 1.0;
             uv += aa;
             uv *= sensorSize / (2.0 * focalLength);

        vec3 viewDirection = normalize(vec3(uv.x, uv.y, -1.0));

        // start with camera forward vector
        vec3 angle = gbufferModelViewInverse[2].xyz;
        float tilt = degrees(asin(angle.y));

        if (abs(angle.y) > 0.99) {
            // use camera up vector when looking completely up/down
            angle = -gbufferModelViewInverse[1].xyz * sign(angle.y);
        }

        float yaw = degrees(atan(angle.z, angle.x) + 180.0);

        viewPosition = Rotate(viewPosition, vec3(1.0, 0.0, 0.0), radians(tilt));
        viewDirection = Rotate(viewDirection, vec3(1.0, 0.0, 0.0), radians(tilt));

        viewPosition = Rotate(viewPosition, vec3(0.0, 1.0, 0.0), radians(yaw));
        viewDirection = Rotate(viewDirection, vec3(0.0, 1.0, 0.0), radians(yaw));
    #endif

    #ifdef USER_DEFINED_COEFFICIENTS
        float mieCoefficient = BetaM_Arbitrary(wavelength) / aerosolScatteringAlbedo;
        float ozoneCoefficient = BetaO_Arbitrary(wavelength);
        float rayleighCoefficient = BetaR_Arbitrary(wavelength);
    #else
        float mieCoefficient = BetaM(wavelength);
        #ifdef PREETHAM_OZONE
            float ozoneCoefficient = PreethamBetaO(wavelength);
        #else
            float ozoneCoefficient = BetaO(wavelength);
        #endif
        float rayleighCoefficient = BetaR(wavelength);
    #endif
    float cloudCoefficient = cloudAbsorption[int(wavelength - 390)];
    #ifndef ENABLE_CLOUDS
        cloudCoefficient = 0.0;
    #endif
    float solarIrradiance = solarIrradiance[int(wavelength - 390.0)];
    vec4 baseAttenuationCoefficients = vec4(rayleighCoefficient, mieCoefficient, ozoneCoefficient, cloudCoefficient);
    #ifndef USE_RAYMARCHED_ATMOSPHERE
        AttenuationCoefficients coefficients;
        coefficients.rayleigh = baseAttenuationCoefficients.x;
        coefficients.aerosol = baseAttenuationCoefficients.y;
        coefficients.ozone = baseAttenuationCoefficients.z;
        coefficients.cloud = baseAttenuationCoefficients.w;
        coefficients.mist = 0.0;
        float atmosphere = PathtraceAtmosphereScattering(viewPosition, viewDirection, sunVector, coefficients, solarIrradiance, wavelength);
    #else
        float atmosphere = RaymarchAtmosphereScattering(viewPosition, viewDirection, sunVector, baseAttenuationCoefficients, solarIrradiance, wavelength);
    #endif
    #if CIE_VERSION == 0
        vec3 simulated = SpectrumToXYZExact_CIE2012(atmosphere / wavelengthPDF, wavelength) * xyzToRGBMatrix_D65;
    #elif CIE_VERSION == 1
        vec3 simulated = SpectrumToXYZExact_CIE1931(atmosphere / wavelengthPDF, wavelength) * xyzToRGBMatrix_D65;
    #endif

    if(any(isnan(simulated))) {
        simulated = vec3(0.0);
    }
    if(any(isinf(simulated))) {
        simulated = vec3(3.4e38);
    }

    vec3 previousColor = texture(colortex0, textureCoordinate).rgb;
    simulationOutput.rgb = mix(previousColor, simulated, 1.0 / (++samples));
    simulationOutput.a   = samples;
}