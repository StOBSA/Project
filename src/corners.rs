use indexmap::IndexSet;

#[derive(Clone)]
pub struct Corners {
    pub included : IndexSet<usize>
}

impl Corners {
    pub fn new() -> Self { Self { included:IndexSet::new() } }

    pub fn iter(&self) -> impl Iterator<Item=usize> + Clone + '_ {
        self.included.iter().map(|&i|i)
    }

    pub fn insert(&mut self, n : usize) {
        self.included.insert(n);
    }

    pub fn remove(&mut self, n : &usize) {
        self.included.remove(n);
    }

    pub fn contains(&self, n : &usize) -> bool{
        self.included.contains(n)
    }
}

impl FromIterator<usize> for Corners {
    fn from_iter<T: IntoIterator<Item = usize>>(iter: T) -> Self {
        let mut corners = Corners::new();
        for i in iter {
            corners.insert(i);
        }
        corners
    }
}

use itertools::Itertools;
impl std::fmt::Debug for Corners {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            f.write_str(format!("{:?}", self.included.iter().collect_vec()).as_str())
        }
    }
// #[derive(Clone)]
// pub struct BinaryCorners {
//     pub included : rug::Integer
// }

// impl BinaryCorners {
//     pub fn new() -> BinaryCorners {
//         BinaryCorners{included:rug::Integer::ZERO}
//     }
    
//     pub fn iter(&self) -> impl Iterator<Item=usize> + Clone + '_ {
//         use rug::Complete;
//         BinaryCornersIterator{
//             corners:self.clone(),
//             index:0,
//             mask :rug::Integer::parse("1").unwrap().complete()
//         }
//     }

//     pub fn insert(&mut self, n : usize) {
//         self.included.set_bit(n as u32, true);
//     }

//     pub fn remove(&mut self, n : &usize) {
//         self.included.set_bit(*n as u32, false);
//     }

//     pub fn contains(&self, n : &usize) -> bool{
//         self.included.get_bit(*n as u32)
//     }
// }

// #[derive(Clone)]
// struct BinaryCornersIterator {
//     mask : rug::Integer, 
//     index : usize,
//     corners : BinaryCorners
// }

// impl Iterator for BinaryCornersIterator {
//     type Item = usize;

//     fn next(&mut self) -> Option<Self::Item> {
//         loop {
//             let result = self.corners.included.get_bit(self.index as u32);
//             self.index += 1;
//             self.mask <<= 1;
//             if result {
//                 return Some(self.index-1)
//             } else if self.mask >= self.corners.included {
//                 return None
//             }
//         }
//     }
// }

// impl FromIterator<usize> for BinaryCorners {
//     fn from_iter<T: IntoIterator<Item = usize>>(iter: T) -> Self {
//         let mut corners = BinaryCorners::new();
//         for i in iter {
//             corners.insert(i);
//         }
//         corners
//     }
// }

// impl std::fmt::Debug for BinaryCorners {
//     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
//         f.write_str(format!("{:?}", self.iter().collect_vec()).as_str())
//     }
// }